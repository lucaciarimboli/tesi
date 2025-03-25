n = 128;                // number of discretization points

For i In {0:n}
  t = 2*Pi*i/n; // Linearly spaced t values from 0 to 2*pi
  Point(i+1) = {
	Sqrt(1 + Cos(t)),
        Sin(t) / 2,
        0
	};
EndFor

// Connect points into a closed boundary
Spline(1) = {1:n+1};
Line(2) = {n+1,1};
Line Loop(3) = {1,2};
Physical Curve("Boundary",1) = {3};   // Assign bounday

// Create a plane surface
Plane Surface(4) = {3};
Physical Surface("Omega",2) = {4};   // Assign domain

Mesh.Algorithm = 1;
Mesh.CharacteristicLengthMax = 0.1;
Mesh.Smoothing = 10;
