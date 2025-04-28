// Parameters:
R0 = 5/3;         // major radius
a = .5;          // minor radius
tau = 0;     // triangularity
eps = a / R0;     // aspect ratio
b = 1 - eps;     // elongation

// --------------------------------------------------------
// GENERATE D-SHAPED BOUNDARY (FOR THE VESSEL)

n = 64;		// number of discretization points

For i In {0:n-1}
	x = -1 + 2*i/n; // Linearly spaced t values from -1 to 1 
	Point(i+1) = {x,
		b/a * Sqrt( (1 - (x - eps/2 * (1 - x^2))^2) / ((1 - eps^2 / 4) * (1 + eps*x)^2 + tau*x * (1 + eps/2 * x)) ),
		0}; // Upper boundary
  
	x = x + 2/n;		// Allows to define once the points {-1,0,0} and {1,0,0}
	Point(n+i+1) = {x,
		-b/a * Sqrt( (1 - (x - eps/2 * (1 - x^2))^2) / ((1 - eps^2 / 4) * (1 + eps*x)^2 + tau*x * (1 + eps/2 * x)) ),
		0}; // Lower boundary
EndFor

// Connect points into curves
Spline(1) = {1:n}; // Upper part
Spline(2) = {n+1:2*n}; // Lower part

// Close the loop
Line(3) = { n, 2*n };
Line(4) = { n+1, 1};
Line Loop(5) = {1, 3, -2, 4};
Physical Curve("Boundary",1) = {1,2,3,4};

// Create a plane surface
Plane Surface(6) = {5};
Physical Surface("Omega") = {6};

// Characteristic lengths for refinement:
h = 5 / n;

// ---------------------------------------------------------
// DEFINE LIMITER POINTS
//Point(1000) = {0.7, 0.0, 0};
//Point(1001) = {-0.7, 0.0, 0};

// Embedd limiter points on the domain
//Point{1000} In Surface{6};
//Point{1001} In Surface{6};
//Physical Point("Limiter Points") = {1000, 1001};

// --------------------------------------------------------
// DEFINE A CLOSED LIMITER CURVE (CIRCULAR)

Point(2000) = {0.0, 0.0, 0};
Point(2001) = {-0.8, 0.0, 0};
Point(2002) = {0.8, 0.0, 0};
Point(2003) = {0.0, 0.8, 0};
Point(2004) = {0.0, -0.8, 0};

Circle(3000) = {2001, 2000, 2003};
Circle(3001) = {2003, 2000, 2002};
Circle(3002) = {2002, 2000, 2004};
Circle(3003) = {2004, 2000, 2001};

// Define a Physical Curve group
Curve Loop(3004) = {3000, 3001, 3002, 3003};
Curve{3000,3001,3002,3003} In Surface{6};
Physical Curve("LimiterCurve") = {3000, 3001, 3002, 3003};

//---------------------------------------------------------
// MESH OVER THE DEFINED DOMAIN

Mesh.Algorithm = 1;
Mesh.CharacteristicLengthMax = h;
Mesh.Smoothing = 10;

// Generate Mesh:
Mesh 2;

// Save Mesh:
Save StrCat("D_mesh_lim.msh");
