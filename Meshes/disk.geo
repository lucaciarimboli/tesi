// Parameters
r = 1.0;      // Circle radius
lc = 0.1;     // Characteristic length of elements

// Define the center point
Point(1) = {0, 0, 0, lc};

// Define boundary poinits on the circle
Point(2) = {r, 0, 0, lc};
Point(3) = {-r, 0, 0, lc};
Point(4) = {0, r, 0, lc};
Point(5) = {0, -r, 0, lc};

// Create circular arcs
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 3};
Circle(4) = {3, 1, 4};

// Create a loop and surface
Curve Loop(5) = {4, 1, 2, 3};
Plane Surface(1) = {5};

// Create physical objects
Physical Curve("Circle", 6) = {4, 3, 2, 1};
Physical Surface("Disc", 2) = {1};
