n = 32;
For i In {0:n-1}
	t = 2*Pi*i/n; // Linearly spaced t values from 0 to 2*pi
	Point(i+1) = { Sqrt(1 + Cos(t)), Sin(t) / 2, 0};
EndFor

// Connect points into a closed boundary
Curve(1) = {1:n};
Curve(2) = {n,1};
Line Loop(3) = {1,2};
Physical Curve("Boundary",1) = {1,2};   // Assign boundary

// Create a plane surface
Plane Surface(4) = {3};
Physical Surface("Omega",2) = {4};   // Assign domain

// Characteristic lengths for refinement:
h = 5 / n;

Mesh.Algorithm = 1;
Mesh.CharacteristicLengthMax = h;
Mesh.Smoothing = 10;

// Generate Mesh:
Mesh 2;

// Save mesh:
Save StrCat("drop_2.msh");
