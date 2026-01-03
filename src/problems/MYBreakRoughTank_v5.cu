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

#include "MYBreakRoughTank_v5.h"
#include "particledefine.h"
#include "GlobalData.h"
#include "cudasimframework.cu"


#define MK_par 2

MYBreakRoughTank_v5::MYBreakRoughTank_v5(GlobalData *_gdata) : Problem(_gdata)
{
	// use planes in general
	const bool use_planes = get_option("use_planes", false);
	// use a plane for the bottom
	const bool use_bottom_plane = get_option("bottom-plane", use_planes);
	// Add objects to the tank
	const bool use_cyl = get_option("cylinder", false);
	// Density diffusion type
	const DensityDiffusionType RHODIFF = get_option("density-diffusion", DELTA-SPH);

	//const bool use_geometries = get_option("use-geometries", true);

	//if (use_bottom_plane && !use_planes)
	//	throw std::invalid_argument("cannot use bottom plane if not using planes");

	// Size and origin of the simulation domain
	lx = 82.5;
	ly = 2.0;
	lz = 3.5;

	// Data for problem setup
	slope_length = 36.0;
	slope2_length = 10.0;
	h_length = 35.5;
	//height = .63;
	height = 3.3;
	//beta = 4.2364*M_PI/180.0;
	//beta = 2.86241*M_PI/180.0;  // bed slope = atan(height/slope_length).
	beta = 1.432*M_PI/180.0;
        beta_2 = 11.30993*M_PI/180.0; //b1/eta for run-up slope

	// add DEM
	//const string dem_file = get_option("dem", "cobble_surface_with_slope.txt");


	SETUP_FRAMEWORK(
		viscosity<SPSVISC>,
		boundary<DUMMY_BOUNDARY>
	).select_options(
		RHODIFF,use_planes,
//		add_flags<ENABLE_DEM|ENABLE_PLANES>()
		add_flags<ENABLE_PLANES>()
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
	simparams()->tend = 20.0f; //seconds
	//simparams()->densityDiffCoeff = 1.0;

	//WaveGage
	if (get_option("gages", false)) {
		add_gage(1, 0.3);
		add_gage(0.5, 0.3);
	}

	// Physical parameters
	H = 2.3;
	set_gravity(-9.81f);
	//setMaxFall(H);

	float r0 = m_deltap;

	auto water = add_fluid( 1000.0f);
	//add_fluid( 1000.0f);
	set_equation_of_state(0, 7.0f, 50.f);
	set_kinematic_visc(0, 1.0e-6);
	set_artificial_visc(0.2f);

	//Wave paddle definition:  location, start & stop times, stroke and frequency (2 \pi/period)
	//paddle_length = .7f;
	paddle_length = 2.3f;
	//paddle_width = m_size.y - 2*r0;
	paddle_width = ly -.2*r0;
	//paddle_tstart=0.5f;
	paddle_tstart=300.0f;
	paddle_origin = make_double3(0.25f, r0, 0.0f);
	//paddle_tend = 30.0f;//seconds
	paddle_tend = 330.0f;
	// The stroke value is given at free surface level H
	float stroke = 0.2;
	// m_mbamplitude is the maximal angular value for paddle angle
	// Paddle angle is in [-m_mbamplitude, m_mbamplitude]
	paddle_amplitude = atan(stroke/(2.0*(H - paddle_origin.z)));
	cout << "\npaddle_amplitude (radians): " << paddle_amplitude << "\n";
	paddle_omega = 2.0*M_PI/0.8;		// period T = 0.8 s

	// Drawing and saving times

	add_writer(VTKWRITER, .25);  //second argument is saving time in seconds

	// Name of problem used for directory creation
	//m_name = "MYBreakRoughTank_v5";

	//GeometryID dem = addDEM(dem_file);
	//addDEM(dem_file, DEM_FMT_ASCII, use_geometries ? FT_NOFILL : FT_BORDER);
	//addDEM(dem_file, DEM_FMT_ASCII, FT_BORDER);

	// Building the geometry
	//const float br = (simparams()->boundarytype == MK_BOUNDARY ? m_deltap/MK_par : r0);
	const int num_layers = (simparams()->boundarytype > SA_BOUNDARY) ?
                simparams()->get_influence_layers() : 1;
        const double box_thickness = (num_layers - 1)*m_deltap;
        const double3 slope_origin = make_double3(paddle_origin.x + h_length, 0, -box_thickness);
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
	rotate(paddle, 0, 0, 0);
	//disableCollisions(paddle);

	double rot_correction1 = sin(beta)*box_thickness;
	//double rot_correction2 = sin(beta_2)*box_thickness;
	if (!use_bottom_plane) {
		//GeometryID bottom = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
		//		Point(h_length, 0, 0), 0, ly, paddle_length);
		//	Vector(slope_length/cos(beta), 0.0, slope_length*tan(beta)));
		GeometryID bottom = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                                slope_origin + make_double3(rot_correction1, 0, (1-cos(beta))*box_thickness),
                                lx - h_length - rot_correction1, ly, box_thickness);
		rotate(bottom, 0, beta, 0);
	//	disableCollisions(bottom);
	//}
	//if (!use_bottom_plane)  {
        //      addPlane(-sin(beta),0,cos(beta), h_length*sin(beta)) ;  //sloping bottom starting at x=h_length
        //      addPlane(-sin(beta_2),0,cos(beta_2), (h_length+slope_length)*sin(beta_2)-cos(beta_2)*tan(beta)*slope_length);
        }

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
                GeometryID bottom = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                        Point(paddle_origin - make_double3(box_thickness, m_deltap, box_thickness)),
                        h_length + box_thickness + rot_correction1, ly, box_thickness);
                setUnfillRadius(bottom, 0.5*m_deltap);

                const double wall_height = paddle_length + box_thickness + (lz - paddle_length)/3.0;
                // close wall
                GeometryID wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                        //Point(m_origin - make_double3(0, box_thickness, box_thickness)),
                        Point(make_double3(0,0,0) - make_double3(0, box_thickness, box_thickness)),
                        lx + paddle_origin.x, box_thickness, wall_height);

                // far wall
                wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                        //Point(m_origin + make_double3(0, ly, -box_thickness)),
                        Point(make_double3(0,0,0) + make_double3(0, ly, -box_thickness)),
                        lx + paddle_origin.x, box_thickness, wall_height);
		// end wall
                wall = addBox(GT_FIXED_BOUNDARY, FT_BORDER,
                        //Point(m_origin + make_double3(0, ly, -box_thickness)),
                        Point(make_double3(0,0,0) + make_double3(lx,0, -box_thickness)),
                        box_thickness, ly,  wall_height);
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

                // this plane corresponds to the initial paddle position, and is only used to cut out
                // the fluid behind the paddle. it will not be an actual geometry
                const double pcx = cos(paddle_amplitude);
                const double pcz = sin(paddle_amplitude);
                const double pcd = paddle_origin.x*pcx + paddle_origin.z*pcz;
                plane = addPlane(pcx, 0, pcz, -pcd, FT_UNFILL);

                setEraseOperation(plane, ET_ERASE_FLUID);
        }
	GeometryID fluid;
	float z = 0;
	int n = 0;
	while (z < H) {
		z = n*(m_deltap+1e-6) + 1.5*r0;    //z = n*m_deltap + 1.5*r0;
		//float x = paddle_origin.x + (z - paddle_origin.z)*tan(amplitude) + 1.0*r0/cos(amplitude);
		float x = paddle_origin.x +r0;
		//float l = h_length + z/tan(beta) - 1.5*r0/sin(beta) - x;
		//float l = h_length;
		float l;
                if (z <= 0.6f) {
                     l = h_length + z/tan(beta) - 1.5*r0/sin(beta) - x;
                //} else if (z <= 0.8f) {
                //     l = h_length + 0.5f/tan(beta) + (z-0.5f)/tan(beta_2) - 1.5*r0/sin(beta_2) - x;
                } else {
                     l = h_length-10;
                }
		fluid = addRect(GT_FLUID, FT_SOLID, Point(x,  r0, z),
				l, ly-2.0*r0);
		n++;
	 }
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
MYBreakRoughTank_v5::moving_bodies_callback(const uint index, Object* object, const double t0, const double t1,
			const float3& force, const float3& torque, const KinematicData& initial_kdata,
			KinematicData& kdata, double3& dx, EulerParameters& dr)
{

    dx= make_double3(0.0);
    kdata.lvel=make_double3(0.0f, 0.0f, 0.0f);
    if (t1> paddle_tstart && t1 < paddle_tend){
       kdata.avel = make_double3(0.0, paddle_amplitude*paddle_omega*sin(paddle_omega*(t1-paddle_tstart)),0.0);
       EulerParameters dqdt = 0.5*EulerParameters(kdata.avel)*kdata.orientation;
       dr = EulerParameters::Identity() + (t1-t0)*dqdt*kdata.orientation.Inverse();
       dr.Normalize();
	   kdata.orientation = kdata.orientation + (t1 - t0)*dqdt;
	   kdata.orientation.Normalize();
	   }
	else {
	   kdata.avel = make_double3(0.0,0.0,0.0);
	   kdata.orientation = kdata.orientation;
	   dr.Identity();
	}
}

//void MYBreakRoughTank_v5::copy_planes(PlaneList &planes)
//{
//	const double w = m_size.y;
//	const double l = h_length + slope_length + slope2_length;

	//  plane is defined as a x + by +c z + d= 0
////	planes.push_back( implicit_plane(0, 0, 1.0, 0) );   //bottom, where the first three numbers are the normal, and the last is d.
//	planes.push_back( implicit_plane(0, 1.0, 0, 0) );   //wall
//	planes.push_back( implicit_plane(0, -1.0, 0, w) ); //far wall
//	planes.push_back( implicit_plane(1.0, 0, 0, 0) );  //end
//	planes.push_back( implicit_plane(-1.0, 0, 0, l) );  //one end
//	if (use_bottom_plane)  {
//		planes.push_back( implicit_plane(-sin(beta),0,cos(beta), h_length*sin(beta)) );  //sloping bottom starting at x=h_length
//	}
	//if (use_bottom_plane)  {
	//	planes.push_back( implicit_plane(-sin(beta_2),0,cos(beta_2), (h_length+slope_length)*sin(beta_2)-cos(beta_2)*tan(beta)*slope_length) );
	//}
//}

#undef MK_par
