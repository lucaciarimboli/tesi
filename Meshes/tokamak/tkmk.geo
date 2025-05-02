//----------------------------------------//
// 	   ARTIFICIAL BOUNDARY		  //
//----------------------------------------//

// Parameters
L = 5.0; // Side length of the square

// Points
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, L, 0};
Point(4) = {0, L, 0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Line Loop and Plane Surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//----------------------------------------//
//	      VACUUM VESSEL		  //
//----------------------------------------//

// Vessel shape parameters:
R0 = 2.0;	// major radius
w = 2.0;	// vessel width
tau = 0;        // vessel triangularity
h = 2.0;	// vessel height	

n = 60; // number of discretization points

eps = 0.5 * w / R0;  // aspect ratio
b = 1 - eps;   	     // elongation

Point(10) = {R0 - w/2, 0, 0};
Point(n+10) = {R0 + w/2, 0, 0};

For i In {1:n-1}
	x = -1 + 2*i/n;
	Point(i+10) = { x + 1 + R0 - w/2,
		      h * 2*b/w*Sqrt((1-(x-eps/2*(1-x^2))^2)/((1-eps^2/4)*(1+eps*x)^2+tau*x*(1+eps/2*x))),
		      0};
EndFor

Spline(10) = {10:n+10}; // Vessel Wall
// Line(11) = {n+10,10};

Curve{10} In Surface{1};

//----------------------------------------//
//		 LIMITER		  //
//----------------------------------------//

//Point(n+11) = {R0 + 2/3*w/2, 0, 0};
//Point(n+12) = {R0 + 2/3*w/2, 0.2, 0};
//Line(12) = {n+11, n+12};
//Curve{12} In Surface{1};	 

//----------------------------------------//
//	          COILS			  //
//----------------------------------------//

// Coils are squares 0.3x0.3
// Coil 1
Point(100) = {R0, 1.4, 0};
Point(101) = {R0+0.3, 1.4, 0};
Point(102) = {R0+0.3, 1.7, 0};
Point(103) = {R0, 1.7, 0};

Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 100};

Line Loop(100) = {100, 101, 102, 103};
Line{100,101,102,103} In Surface{1};

// Coil 2:
Point(110) = {R0+1.5, 0.5, 0};
Point(111) = {R0+1.8, 0.5, 0};
Point(112) = {R0+1.5, 0.8, 0};
Point(113) = {R0+1.8, 0.8, 0};

Line(110) = {110, 111};
Line(111) = {111, 112};
Line(112) = {112, 113};
Line(113) = {113, 110};

Line Loop(110) = {110, 111, 112, 113};
Line{110,111,112,113} In Surface{1};

//----------------------------------------//
// 	     PHYSICAL QUANTITIES	  //
//----------------------------------------//

Physical Surface("AirRegion") = {1};

Physical Curve("Gamma0") = {2, 3, 4};
Physical Curve("Gamma1") = {1};
Physical Curve("VesselWall") = {10};

// Physical Curve("Limiter") = {12};

Physical Curve("Coil1") = {100,101,102,103};
Physical Curve("Coil2") = {110,111,112,113};

//----------------------------------------//
//	       GENERATE MESH		  //
//----------------------------------------//

Mesh.Algorithm = 1;
// Mesh.CharacteristicLengthMax = 0.5;
Mesh.Smoothing = 10;

Mesh 2; // generate mesh
Save StrCat("tkmk.msh");  // save mesh	
