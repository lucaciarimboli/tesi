h = 0.05;
n = 64;

For i In {0:n-1}
	t = 2*Pi*i/n; // Linearly spaced t values from 0 to 2*pi
	Point(i+1) = {Sqrt(1 + Cos(t)), Sin(t) / 2, 0};
EndFor

// Connect points into a closed boundary
Curve(1) = {1:n-1};
Curve(2) = {n-1,1};
Line Loop(3) = {1,2};
Physical Curve("Boundary",1) = {1,2};   // Assign boundary

// Create a plane surface
Plane Surface(4) = {3};
Physical Surface("Omega",2) = {4};   // Assign domain

Mesh.Algorithm = 4;
Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthMax = 1/n;
Mesh.Smoothing = 10;
Mesh 2;
Save "drop_2.msh";
