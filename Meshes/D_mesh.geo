// Parameters:
R0 = 5/3;         // major radius
a = .5;          // minor radius
lambda = 0;     // triangularity
eps = a / R0;     // aspect ratio
b = 1 - eps;     // elongation

// Function for the boundary curve
yplus(x) = b/a * Sqrt( (1 - (x - eps/2 * (1 - x^2))^2) / ((1 - eps^2 / 4) * (1 + eps*x)^2 + lambda*x * (1 + eps/2 * x)) );
Function yminus(x) = -yplus(x);

// Discretize the boundary
For i In {0:n}
  xi = -1 + 2*i/n; // Linearly spaced x values from -1 to 1
  Point(i+1) = {xi, yplus(xi), 0}; // Upper boundary
  Point(n+i+2) = {xi, yminus(xi), 0}; // Lower boundary
EndFor

// Connect points into boundary curves
Spline(1) = {1:n+1}; // Upper part
Spline(2) = {n+2:2*n+2}; // Lower part

// Close the loop
Line(3) = {1, 2*n+2};
Line Loop(4) = {1, 3, 2};

// Create a plane surface
Plane Surface(5) = {4};

// Mesh settings
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.1;

// Generate the mesh
Mesh 2;

