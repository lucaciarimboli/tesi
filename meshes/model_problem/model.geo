//----------------------------------//
//		REFINMENT PARAMETERS		//
//----------------------------------//
lc1 = 0.3;    // Boundary size
lc2 = 0.1;  // Coils size

//----------------------------------//
//	           BOUNDARY             //
//----------------------------------//
Point(1) = {0.01, -2, 0, lc1};
Point(2) = {3, -2, 0, lc1};
Point(3) = {3, 2, 0, lc1};
Point(4) = {0.01, 2, 0, lc1};

// Point(1) = {0.01, -40, 0, lc1};
// Point(2) = {60, -40, 0, lc1};
// Point(3) = {60, 40, 0, lc1};
// Point(4) = {0.01, 40, 0, lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

//----------------------------------//
//	       	     COILS             	//
//----------------------------------//
// PF COIL 1
Point(5) = {1, 0.5, 0, lc2};
Point(6) = {2, 0.5, 0, lc2};
Point(7) = {2, 1.5, 0, lc2};
Point(8) = {1, 1.5, 0, lc2};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

Line Loop(2) = {5,6,7,8};	

// PF COIL 2
Point(9) = {0.5,-1.5,0,lc2};
Point(10) = {1.5,-1.5,0,lc2};
Point(11) = {1.5, -0.5, 0, lc2};
Point(12) = {0.5, -0.5, 0 , lc2};

Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,9};

Line Loop(3) = {9,10,11,12};

//----------------------------------//
//	        DOMAIN REGIONS	        //
//----------------------------------//
Plane Surface(1) = {1,-2,-3}; // Air
Plane Surface(2) = {2};
Plane Surface(3) = {3};

//----------------------------------//
//	       PHYSICAL GROUPS	        //
//----------------------------------//
Physical Surface("Air",1) = {1};
Physical Surface("Coil 1",2) = {2};
Physical Surface("Coil 2",3) = {3};

//Physical Curve("Boundary",4) = {1,2,3,4};
Physical Curve("Coil 1", 5) = {5,6,7,8};
Physical Curve("Coil 2", 6) = {9,10,11,12};
Physical Curve("Dirichlet", 7) = {4};
Physical Curve("Neumann", 8) = {1,2,3};

//----------------------------------//
//	        GENERATE MESH	        //
//----------------------------------//
Mesh.Algorithm = 5;
Mesh.MshFileVersion = 2.2;
Mesh 2;
Save "model.msh";