n = 128;                // number of discretization points

For i In {0:n-1}
  t = 2*Pi*i/n; // Linearly spaced t values from 0 to 2*pi
  Point(i+1) = {
	Sqrt(1 + Cos(t)),
        Sin(t) / 2,
        0
	};
EndFor

// Connect points into a closed boundary
Spline(1) = {1:n};
Line(2) = {n,1};
Line Loop(3) = {1,2};
Physical Curve("Boundary",1) = {1,2};   // Assign boundary

// Create a plane surface
Plane Surface(4) = {3};
Physical Surface("Omega",2) = {4};   // Assign domain

// Characteristic lengths for refinement:
h_values[] = {0.8, 0.5, 0.3, 0.1};

For i In {0:#h_values[]-1}

	Mesh.Algorithm = 1;
	Mesh.CharacteristicLengthMax = h_values[i];
	Mesh.Smoothing = 10;

	// Generate Mesh:
    	Mesh 2;

    	// Save mesh:
    	Save StrCat("drop_", Sprintf("%g", i+1), ".msh");
EndFor
