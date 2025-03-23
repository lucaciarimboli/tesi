// Parameters
radius = 1.0;  // Circle radius
lc = 0.05;     // Characteristic length of elements

// Define the center point
Point(1) = {0, 0, 0, lc};

// Define boundary points on the circle
Point(2) = {radius, 0, 0, lc};
Point(3) = {0, radius, 0, lc};
Point(4) = {-radius, 0, 0, lc};
Point(5) = {0, -radius, 0, lc};

// Create circular arcs
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Create a loop and surface
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Mesh generation
Mesh 2;  // Generate a 2D mesh

// Save the mesh
Save "circle_mesh.msh";

