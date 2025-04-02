// Parameters:
R0 = 5/3;         // major radius
a = .5;          // minor radius
tau = 0;     // triangularity
eps = a / R0;     // aspect ratio
b = 1 - eps;     // elongation

n = 128;		// number of discretization points

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
Physical Surface("Omega",2) = {6};

// Characteristic lengths for refinement:
h_values[] = {0.8, 0.5, 0.3, 0.1};

For i In {0:#h_values[]-1}

        Mesh.Algorithm = 1;
        Mesh.CharacteristicLengthMax = h_values[i];
        Mesh.Smoothing = 10;

        // Generate Mesh:
        Mesh 2;

        // Save mesh:
        Save StrCat("D_mesh_", Sprintf("%g", i+1), ".msh");
EndFor

