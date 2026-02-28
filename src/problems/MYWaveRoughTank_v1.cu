/*  Copyright (c) 2011-2019 INGV, EDF, UniCT, JHU

    Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Catania, Italy
    Électricité de France, Paris, France
    Università di Catania, Catania, Italy
    Johns Hopkins University, Baltimore (MD), USA

    This file is part of GPUSPH. Project founders:
        Alexis Hérault, Giuseppe Bilotta, Robert A. Dalrymple,
        Eugenio Rustico, Ciro Del Negro
    For a full list of authors and project partners, consult the logs
    and the project website <https://www.gpusph.org>

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdexcept>

#include "MYWaveRoughTank_v1.h"
#include "particledefine.h"
#include "GlobalData.h"
#include "cudasimframework.cu"


#define MK_par 2

MYWaveRoughTank_v1::MYWaveRoughTank_v1(GlobalData *_gdata) : Problem(_gdata)
{
	// use planes in general
	const bool use_planes = get_option("use_planes", false);
	// use a plane for the bottom
	const bool use_bottom_plane = get_option("bottom-plane", use_planes);
	// Add objects to the tank
	const bool use_cyl = get_option("cylinder", false);
	// Density diffusion type
	const DensityDiffusionType RHODIFF = get_option("density-diffusion", DELTA_SPH);

	const bool use_geometries = get_option("use-geometries", true);

	//if (use_bottom_plane && !use_planes)
	//	throw std::invalid_argument("cannot use bottom plane if not using planes");

	// Size and origin of the simulation domain
	lx = 86.5;
	ly = 2.0;
	lz = 4.5;

	// Data for problem setup
	slope_length = 56.0;
	slope2_length = 15.0;
	h_length = 15.5;
	//height = .63;
	height = 3.8;
	//beta = 4.2364*M_PI/180.0;
	//beta = 2.86241*M_PI/180.0;  // bed slope = atan(height/slope_length).
	beta = 1.432*M_PI/180.0;
        beta_2 = 11.30993*M_PI/180.0; //b1/eta for run-up slope

	// add DEM
	//const string dem_file = get_option("dem", "cobble_surface_with_slope_v3.txt");
	//const string dem_file = get_option("dem", "cobble_surface_with_slope_v6_res001_2.0fac.txt");
	//const string dem_file = get_option("dem", "cobble_surface_with_slope_v4_0.3fac.txt");
	const string dem_file = get_option("dem", "preprocessing/cobble_slope_design/cobble_surface_with_slope_wave_v1_fac1.0.txt");


	SETUP_FRAMEWORK(
        //        rheology<NEWTONIAN>,
        //        turbulence_model<ARTIFICIAL>,
		viscosity<SPSVISC>,
		boundary<DUMMY_BOUNDARY>
	).select_options(
		RHODIFF,use_geometries,
		add_flags<ENABLE_DEM|ENABLE_PLANES>()
//		add_flags<ENABLE_PLANES>()
		//add_flags<ENABLE_DEM | ENABLE_PLANES>
	);

	// Allow user to set the MLS frequency at runtime. Default to 0 if density
	// diffusion is enabled, 10 otherwise
	const int mlsIters = get_option("mls",
		(simparams()->densitydiffusiontype != DENSITY_DIFFUSION_NONE) ? 0 : 10);

	if (mlsIters > 0)
		addFilter(MLS_FILTER, mlsIters);


	//m_size = make_double3(lx, ly, lz);
	//m_origin = make_double3(0, 0, 0);
	//if (use_cyl) {
	//	m_origin.z -= 2.0*height;
	//	m_size.z += 2.0*height;
	//}

	//addFilter(SHEPARD_FILTER, 20); // or MLS_FILTER

	if (get_option("testpoints", false)) {
		addPostProcess(TESTPOINTS);
	}

	//// use a plane for the bottom
	//use_bottom_plane = 1;  //1 for plane; 0 for particles

	// SPH parameters
	set_deltap(0.03f);  //0.005f;
	//set_timestep(0.0001);
	simparams()->dtadaptfactor = 0.2;
	simparams()->buildneibsfreq = 10;
	simparams()->tend = 30.0f; //seconds
	simparams()->densityDiffCoeff = 1.0;

	//WaveGage
	if (get_option("gages", false)) {
		add_gage(1, 0.3);
		add_gage(0.5, 0.3);
	}

	// Physical parameters
	H = 1.5;
	float water_height = 0.8;
	set_gravity(-9.81f);
	//setMaxFall(H);
	//setMaxParticleSpeed(7.0);

	float r0 = m_deltap;

	auto water = add_fluid( 1000.0f);
	//add_fluid( 1000.0f);
	set_equation_of_state(0, 7.0f, 50.f);
	set_kinematic_visc(0, 1.0e-6);
	set_artificial_visc(0.2f);

	//Wave paddle definition:  location, start & stop times, stroke and frequency (2 \pi/period)
	//paddle_length = .7f;
	paddle_length = 4.3f;
	//paddle_width = m_size.y - 2*r0;
	paddle_width = ly - 10*r0;
	paddle_tstart=0.5f;
	paddle_origin = make_double3(5*r0, 6*r0, 4*r0);
	paddle_tend = 30.0f;
	// The stroke value is given at free surface level H
	// float stroke = 0.2;
	// m_mbamplitude is the maximal angular value for paddle angle
	// Paddle angle is in [-m_mbamplitude, m_mbamplitude]
	//paddle_amplitude = atan(stroke/(2.0*(H - paddle_origin.z)));
	paddle_amplitude = 2;
	cout << "\npaddle_amplitude (radians): " << paddle_amplitude << "\n";
	paddle_omega = 2.0*M_PI/10.0f;		// period T = 0.8 s
	paddle_period = 10.0f;

	// Drawing and saving times

	add_writer(VTKWRITER, .25);  //second argument is saving time in seconds

	// Name of problem used for directory creation
	//m_name = "MYWaveRoughTank_v1";

	//GeometryID dem = addDEM(dem_file);
	//addDEM(dem_file, DEM_FMT_ASCII, use_geometries ? FT_NOFILL : FT_BORDER);
	addDEM(dem_file, DEM_FMT_ASCII, FT_BORDER);

	// Building the geometry
	//const float br = (simparams()->boundarytype == MK_BOUNDARY ? m_deltap/MK_par : r0);
	//const int num_layers = (simparams()->boundarytype > SA_BOUNDARY) ?
        //        simparams()->get_influence_layers() : 1;
	const int num_layers = 5;
        const double box_thickness = (num_layers - 1)*m_deltap;
        const double3 slope_origin = make_double3(paddle_origin.x + h_length, 0, -box_thickness);
	const double3 slope_origin_2 = make_double3(paddle_origin.x + h_length + slope_length, 0, slope_length*tan(beta)-box_thickness);
        setDynamicBoundariesLayers(num_layers);

	setPositioning(PP_CORNER);

	//GeometryID experiment_box = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
	//Point(0, 0, 0), h_length + slope_length + slope2_length,ly, height);
	//disableCollisions(experiment_box);

	const float amplitude = -paddle_amplitude ;
	GeometryID paddle = addBox(GT_MOVING_BODY, FT_BORDER,
		Point(paddle_origin- make_double3(box_thickness, 0, 0)),
		box_thickness, paddle_width, paddle_length);
	//rotate(paddle, 0,-amplitude, 0);
	//rotate(paddle, 0, 0, 0);
	//disableCollisions(paddle);

	double rot_correction1 = sin(beta)*box_thickness;
	double rot_correction2 = sin(beta_2)*box_thickness;
	//if (!use_bottom_plane) {
		//GeometryID bottom = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
		//		Point(h_length, 0, 0), 0, ly, paddle_length);
		//	Vector(slope_length/cos(beta), 0.0, slope_length*tan(beta)));
	//	GeometryID bottom = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
        //                        slope_origin + make_double3(rot_correction1, 0, (1-cos(beta))*box_thickness),
        //                        (slope_length - rot_correction1)/cos(beta), ly, box_thickness);
	//	rotate(bottom, 0, beta, 0);
	//	GeometryID bottom2 = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
        //                        slope_origin_2 + make_double3(rot_correction2, 0, (1-cos(beta_2))*box_thickness),
        //                        (slope2_length - rot_correction2)/cos(beta_2), ly, box_thickness);
	//	rotate(bottom2, 0, beta_2, 0);
	//	disableCollisions(bottom);
	//}
	//if (!use_bottom_plane)  {
        //      addPlane(-sin(beta),0,cos(beta), h_length*sin(beta)) ;  //sloping bottom starting at x=h_length
        //      addPlane(-sin(beta_2),0,cos(beta_2), (h_length+slope_length)*sin(beta_2)-cos(beta_2)*tan(beta)*slope_length);
        //}

	//GeometryID dem = addDEM(dem_file);
	if (use_planes) {
                const double w = m_size.y;
                //const double l = h_length + slope_length;
		const double l = h_length + slope_length + slope2_length;

                addPlane(0, 0, 1, 0);  //bottom, where the first three numbers are the normal, and the last is d.
                addPlane(0, 1, 0, 0);  //wall
                addPlane(0, -1, 0, w); //far wall
                addPlane(1.0, 0, 0, 0);   //end
                addPlane(-1.0, 0, 0, l);  //one end
	} else {
                // flat bottom rectangle (before the slope begins)
        //        GeometryID bottom = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
        //                Point(paddle_origin - make_double3(box_thickness, m_deltap, box_thickness)),
        //                h_length + box_thickness + rot_correction1, ly, box_thickness);
        //        setUnfillRadius(bottom, 0.5*m_deltap);

        //        const double wall_height = paddle_length + box_thickness + (lz - paddle_length)/3.0;
	//	cout << "\nwall height: " << wall_height << "\n";
                // close wall
                GeometryID wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                //        //Point(m_origin - make_double3(0, box_thickness, box_thickness)),
                        Point(make_double3(0,0,0) - make_double3(box_thickness+m_deltap, -m_deltap, -m_deltap)),
                        lx + 2*box_thickness, box_thickness, lz-water_height);

                // far wall
                wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                //        //Point(m_origin + make_double3(0, ly, -box_thickness)),
                        Point(make_double3(0,0,0) + make_double3(0, ly-3*m_deltap-0.02, -m_deltap)),
                        lx + 2*box_thickness, box_thickness, lz);
		// end (right) wall
                wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                //        //Point(m_origin + make_double3(0, ly, -box_thickness)),
                        Point(make_double3(0,0,0) + make_double3(h_length+slope_length+slope2_length-m_deltap,0, -box_thickness)),
        		box_thickness, ly,  lz);
		// left wall
		//wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                //        //Point(m_origin + make_double3(0, ly, -box_thickness)),
                //        Point(make_double3(0,0,0) + make_double3(-box_thickness-m_deltap,0, water_height+m_deltap)),
                //        box_thickness, ly,  lz-water_height);
		//top cover
		wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                        //Point(m_origin + make_double3(0, ly, -box_thickness)),
                        Point(make_double3(0,0,0) + make_double3(-box_thickness, -m_deltap, lz)),
                        lx + 2*box_thickness, ly + 2*m_deltap, box_thickness);
        }
	GeometryID fluid;
        float z = 0;
        int n = 0;
        while (z < H) {
                z = n*(m_deltap+1e-6) + 5*r0;    //z = n*m_deltap + 1.5*r0;
                //z = n*(m_deltap+1e-6) + water_height;
                //float x = paddle_origin.x + (z - paddle_origin.z)*tan(amplitude) + 1.0*r0/cos(amplitude);
                //float x = paddle_origin.x +r0;
                float x = 0.0;
                //float l = h_length + z/tan(beta) - 1.5*r0/sin(beta) - x;
                //float l = h_length;
                float l;
                //if (z <= 0.6f) {
                     l = h_length + z/tan(beta) - 5*r0/sin(beta) - x;
                //} else if (z <= 0.8f) {
                //     l = h_length + 0.5f/tan(beta) + (z-0.5f)/tan(beta_2) - 1.5*r0/sin(beta_2) - x;
                //} else {
                //     l = h_length-10;
                //}
                fluid = addRect(GT_FLUID, FT_SOLID, Point(x+6*r0,  6*r0, z),
                                l, ly-11*r0);
                n++;
         }
	// these planes are used at least for cutting, so they are always defined
        {
                // sloping bottom as a plane. if use_bottom_plane, then it will be
                // an actual geometry; if !use_bottom_plane, it will only be used
                // to unfill the fluid (since the sloping box would not be sufficient
                // to remove all of the fluid below)
                GeometryID plane = addPlane(-sin(beta), 0, cos(beta), slope_origin.x*sin(beta),
                        use_bottom_plane ? FT_NOFILL : FT_UNFILL);

                setEraseOperation(plane, ET_ERASE_FLUID);

                // this plane cuts the lateral walls below the sloping ground
                plane = addPlane(-sin(beta), 0, cos(beta),
                        slope_origin.x*sin(beta) + 2*(m_deltap + box_thickness*cos(beta)),
                        FT_UNFILL);

                setEraseOperation(plane, ET_ERASE_BOUNDARY);

		plane = addPlane(-sin(beta_2), 0, cos(beta_2), slope_origin_2.x*sin(beta_2),
                        use_bottom_plane ? FT_NOFILL : FT_UNFILL);

                setEraseOperation(plane, ET_ERASE_FLUID);

                // this plane cuts the lateral walls below the sloping ground
                plane = addPlane(-sin(beta_2), 0, cos(beta_2),
                        slope_origin_2.x*sin(beta_2) - slope_length*sin(beta) + 2*(m_deltap + box_thickness*cos(beta_2)),
                        FT_UNFILL);

                setEraseOperation(plane, ET_ERASE_BOUNDARY);

                // this plane corresponds to the initial paddle position, and is only used to cut out
                // the fluid behind the paddle. it will not be an actual geometry
                const double pcx = cos(paddle_amplitude);
                const double pcz = sin(paddle_amplitude);
                const double pcd = paddle_origin.x*pcx + paddle_origin.z*pcz;
                //plane = addPlane(pcx, 0, pcz, -pcd, FT_UNFILL);

                //setEraseOperation(plane, ET_ERASE_FLUID);
        }
	//float water_height = 0.8;
	//addDEMFluidBox(water_height);
	//addExtraWorldMargin(5*m_deltap);
	//GeometryID fluid;
	//  moved upper sections
	//}
// activate the solid obstacle
	//const uint NUM_OBSTACLES = 100;
	//const double Y_DISTANCE = ly / (NUM_OBSTACLES + 1);
        // // rotation angle
        //const double Z_ANGLE = M_PI / 4;
        //for (uint i = 0; i < NUM_OBSTACLES; i++) {
        //        // Obstacle is of type GT_MOVING_BODY, although the callback is not even implemented, to
        //        // make the forces feedback available
	//	const double dx = i*0.1;
	//	const double height_obstacle = (obstacle_xpos+dx-h_length-slope_length+obstacle_side)*tan(beta_2)+slope_length*tan(beta);
	//	//cout << "\nheight_obstacle: " << height_obstacle << "\n";

        //        //GeometryID obstacle = addBox(GT_MOVING_BODY, FT_BORDER,
	//	GeometryID obstacle = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
        //                Point(obstacle_xpos+dx, r0, 0),
        //                        obstacle_side, ly-2.0*r0, height_obstacle );
        //        if (ROTATE_OBSTACLE) {
        //                rotate(obstacle, 0, 0, Z_ANGLE);
        //                // until we'll fix it, the rotation centers are always the corners
        //                // shift(obstacle, 0, obstacle_side/2, 0);
        //        }
        //        // enable force feedback to measure forces
        //        //enableFeedback(obstacle);
	//	disableCollisions(obstacle);
        //}
// add DEM
	//GeometryID dem = addDEM(dem_file);

	if (hasPostProcess(TESTPOINTS)) {
		Point pos = Point(0.5748, 0.1799, 0.2564, 0.0);
		addTestPoint(pos);
		pos = Point(0.5748, 0.2799, 0.2564, 0.0);
		addTestPoint(pos);
		pos = Point(1.5748, 0.2799, 0.2564, 0.0);
		addTestPoint(pos);
	}

	if (use_cyl) {
		setPositioning(PP_BOTTOM_CENTER);
		Point p[10];
		p[0] = Point(h_length + slope_length/(cos(beta)*10), ly/2., 0);
		p[1] = Point(h_length + slope_length/(cos(beta)*10), ly/6.,  0);
		p[2] = Point(h_length + slope_length/(cos(beta)*10), 5*ly/6, 0);
		p[3] = Point(h_length + slope_length/(cos(beta)*5), 0, 0);
		p[4] = Point(h_length + slope_length/(cos(beta)*5), ly/3, 0);
		p[5] = Point(h_length + slope_length/(cos(beta)*5), 2*ly/3, 0);
		p[6] = Point(h_length + slope_length/(cos(beta)*5), ly, 0);
		p[7] = Point(h_length + 3*slope_length/(cos(beta)*10), ly/6, 0);
		p[8] = Point(h_length + 3*slope_length/(cos(beta)*10), ly/2, 0);
		p[9] = Point(h_length+ 3*slope_length/(cos(beta)*10), 5*ly/6, 0);
		p[10] = Point(h_length+ 4*slope_length/(cos(beta)*10), ly/2, 0);

		for (int i = 0; i < 11; i++) {
			GeometryID cyl = addCylinder(GT_FIXED_BOUNDARY, FT_BORDER,
				p[i], .025, height);
			disableCollisions(cyl);
			setEraseOperation(cyl, ET_ERASE_FLUID);
		}
	}
}


void
MYWaveRoughTank_v1::moving_bodies_callback(const uint index, Object* object, const double t0, const double t1,
			const float3& force, const float3& torque, const KinematicData& initial_kdata,
			KinematicData& kdata, double3& dx, EulerParameters& dr)
{
    dx= make_double3(0.0);
    kdata.lvel=make_double3(0.0f, 0.0f, 0.0f);
    cout << "\nmove.lvel.x: " << t1 << "\n";
    if (t1> paddle_tstart && t1 < paddle_tend){
       //kdata.avel = make_double3(0.0, paddle_amplitude*paddle_omega*sin(paddle_omega*(t1-paddle_tstart)),0.0);
	//kdata.lvel = make_double3(paddle_amplitude*2.0*M_PI*sin(2.0*M_PI*(t1-paddle_tstart)/paddle_period),0,0);
	float lvel0=0;
	float paddle_shift=7.5; //seconds
    	float wt = paddle_omega * (t1-paddle_tstart+paddle_shift);
	float a1 = -3.807040403e-16;
        float b1 = 1.833777260e+00;
        float a2 = -2.667460375e-01;
        float b2 = -1.109393431e-16;
        float a3 = -8.104411239e-17;
        float b3 = -5.084284076e-02;
        float a4 = 1.080287072e-02;
        float b4 = 3.555081983e-17;
        float a5 = -4.452772794e-18;
        float b5 = 2.479574682e-03;
        float a6 = -6.217843014e-04;
        float b6 = 1.429394875e-17;
        float a7 = 5.104491855e-17;
        float b7 = -1.709078608e-04;
        float a8 = 4.979581269e-05;
        float b8 = -2.549119733e-18;
	lvel0 += paddle_amplitude * (1*paddle_omega) * (-a1*sinf(1*wt) + b1*cosf(1*wt));
	lvel0 += paddle_amplitude * (2*paddle_omega) * (-a2*sinf(2*wt) + b2*cosf(2*wt));
	lvel0 += paddle_amplitude * (3*paddle_omega) * (-a3*sinf(3*wt) + b3*cosf(3*wt));
	lvel0 += paddle_amplitude * (4*paddle_omega) * (-a4*sinf(4*wt) + b4*cosf(4*wt));
        lvel0 += paddle_amplitude * (5*paddle_omega) * (-a5*sinf(5*wt) + b5*cosf(5*wt));
        lvel0 += paddle_amplitude * (6*paddle_omega) * (-a6*sinf(6*wt) + b6*cosf(6*wt));
        lvel0 += paddle_amplitude * (7*paddle_omega) * (-a7*sinf(7*wt) + b7*cosf(7*wt));
        lvel0 += paddle_amplitude * (8*paddle_omega) * (-a8*sinf(8*wt) + b8*cosf(8*wt));
	kdata.lvel = make_double3(lvel0, 0, 0);
       //EulerParameters dqdt = 0.5*EulerParameters(kdata.avel)*kdata.orientation;
       //dr = EulerParameters::Identity() + (t1-t0)*dqdt*kdata.orientation.Inverse();
       //dr.Normalize();
	//   kdata.orientation = kdata.orientation + (t1 - t0)*dqdt;
	//   kdata.orientation.Normalize();
	dx.x = (t1 - t0) * kdata.lvel.x;
	cout << "\nkdata.lvel.x: " << kdata.lvel.x << "\n";
    }
    else {
	   //kdata.avel = make_double3(0.0,0.0,0.0);
	   //kdata.orientation = kdata.orientation;
	   //dr.Identity();
	   dx.x = 0;
    }
}

#undef MK_par
