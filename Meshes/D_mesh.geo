// Parameters:
R0 = 5/3;         // major radius
a = .5;          // minor radius
tau = 0;     // triangularity
eps = a / R0;     // aspect ratio
b = 1 - eps;     // elongation

n = 128;		// number of discretization points

For i In {0:n}
  x = -1 + 2*i/n; // Linearly spaced t values from -1 to 1 
  Point(i+1) = {x,
	b/a * Sqrt( (1 - (x - eps/2 * (1 - x^2))^2) / ((1 - eps^2 / 4) * (1 + eps*x)^2 + tau*x * (1 + eps/2 * x)) ),
	0}; // Upper boundary
  Point(n+i+2) = {x,
	-b/a * Sqrt( (1 - (x - eps/2 * (1 - x^2))^2) / ((1 - eps^2 / 4) * (1 + eps*x)^2 + tau*x * (1 + eps/2 * x)) ),
	0}; // Lower boundary
EndFor

// Connect points into curves
Spline(1) = {1:n+1}; // Upper part
Spline(2) = {n+2:2*n+2}; // Lower part

// Close the loop
Line(3) = {n+2,1};
Line(4) = {n+1, 2*n+2};
Line Loop(5) = {1, 4, -2, 3};
Physical Line("Boundary") = {5};

// Create a plane surface
Plane Surface(6) = {5};
Physical Surface("Omega") = {6};

Mesh.Algorithm = 1;
//Mesh.Smoothing = 100;

