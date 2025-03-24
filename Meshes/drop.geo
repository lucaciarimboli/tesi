n = 128;                // number of discretization points

For i In {0:n}
  t = 2*Pi*i/n; // Linearly spaced t values from 0 to 2*pi
  Point(i+1) = {
	Sqrt(1 + Cos(t)),
        Sin(t) / 2,
        0
	};
EndFor

// Connect points into curves
Spline(1) = {1:n+1};
Line(2) = {n+1,1};
Line Loop(3) = {1, 2};
Physical Line("Boundary") = {3};

// Create a plane surface
Plane Surface(4) = {3};
Recombine Surface {4};
Physical Surface("Omega") = {4};
