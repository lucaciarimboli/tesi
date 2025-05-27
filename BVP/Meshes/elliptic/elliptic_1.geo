n = 16;  // Number of discretization points

// Define R limits
R_min = 2/3 * Sqrt(3*Pi);
R_max = 2/3 * Sqrt(6*Pi);
delta_R = (R_max - R_min) / (n - 1);  // Step size

For i In {0:n-1}
  R = R_min + i * delta_R;
  cos_val = Cos(R^2 / 2);
  
  // Ensure the argument of arccos is in [-1,1]
  If (2*cos_val == 0)
    Printf("Skipping R = %g due to division by zero.", R);
  ElseIf (Abs(1 / (2 * cos_val)) > 1)
    Printf("Skipping R = %g due to invalid arccos argument.", R);
  Else
    Z1 = Acos( 1 / ( 2 * cos_val ) );
    Z2 = 2*Pi - Z1;

    Point(i+1) = {R, Z1, 0}; // Upper boundary
    Point(n+i+1) = {R, Z2, 0}; // Lower boundary
  EndIf
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
h = 5 / n;

// Mesh settings
Mesh.Algorithm = 1;
Mesh.CharacteristicLengthMax = h;
Mesh.Smoothing = 10;

// Generate Mesh:
Mesh 2;

// Save mesh:
Save StrCat("elliptic_1.msh");
