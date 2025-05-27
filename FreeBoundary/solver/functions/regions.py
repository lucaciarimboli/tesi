from firedrake import *

def initialize_plasma_mask(Vessel, V, x, y):
    x0 = Vessel[0]
    y0 = Vessel[1]
    r = Vessel[2]

    plasma_mask = Function(V)
    plasma_mask.interpolate(0.5 + 0.5 * tanh((r**2 - (x - x0)**2 - (y - y0)**2) / 0.01))

    return plasma_mask

def define_vessel_mask(Vessel, j_cv, V, x, y):
    x0 = Vessel[0]
    y0 = Vessel[1]
    r = Vessel[2]
    thickness = Vessel[3]

    vessel_mask = Function(V)
    vessel_mask.interpolate(conditional((x - x0)**2 + (y - y0)**2 < (r + thickness)**2, 1.0, 0.0))
    vessel_mask.assign(Constant(j_cv) * (vessel_mask - 0.5))

    return vessel_mask
    
def define_coils_mask(I, coils, V, x, y):
    n = len(I)

    if n != len(coils):
        raise ValueError("Number of coils and currents do not match!")
    
    coils_mask = Function(V)

    for i in range(n):
        [x_min, x_max, y_min, y_max] = coils[i]

        S = (x_max - x_min) * (y_max - y_min)

        coils_mask.interpolate(conditional(
            And( And(x_min <= x, x <= x_max), And(y_min <= y, y <= y_max) ),
            I[i] / S,
            0.0
        ))

    return coils_mask

def update_plasma_mask(psi, limiter_points, psi_mask):
    # Find the maximum value of psi at the limiter points:
    psi_vals = []
    for (x_pt, y_pt) in limiter_points:
        try:
            psi_vals.append(psi.at((x_pt, y_pt)))
        except:
            raise ValueError(f"Point {(x_pt, y_pt)} outside domain!")
    psi0 = max(psi_vals)

    # Update the plasma boundary mask:
    epsilon = 0.01
    psi_mask.interpolate(0.5 + 0.5 * tanh((psi - psi0) / (epsilon * psi0)))

    return psi0