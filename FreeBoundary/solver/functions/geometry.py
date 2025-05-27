import gmsh

def generate_mesh(params: dict, msh_path: str):
    """
    Use Gmsh Python API to build and save a Firedrake-compatible .msh mesh.
    params: dict with keys:
        - x0, y0: vessel center
        - R: vessel inner radius
        - thickness: vessel wall thickness
        - coil_positions: list of (x_min, x_max, y_min, y_max)
        - domain_size: (width, height) of the outer box (default (1.0, 1.0))
        - mesh_size_min: minimum mesh size (default 0.01)
        - mesh_size_max: maximum mesh size (default 0.05)
        - vacuum_mesh_size: mesh size in vacuum region (default 0.02)
    msh_path: output .msh file path
    """
    gmsh.initialize()
    gmsh.model.add("tokamak")

    # Vessel parameters
    x0 = params["x0"]
    y0 = params["y0"]
    r_inner = params["R"]
    thickness = params["thickness"]
    r_outer = r_inner + thickness

    # Coils
    coil_positions = params["coils"]

    # Domain size
    domain_w, domain_h = params.get("domain_size", (1.0, 1.0))

    # Mesh options
    mesh_size_min = params.get("mesh_size_min", 0.01)
    mesh_size_max = params.get("mesh_size_max", 0.05)
    vacuum_mesh_size = params.get("vacuum_mesh_size", 0.02)

    # Vessel geometry
    outer_disk = gmsh.model.occ.addDisk(x0, y0, 0, r_outer, r_outer)
    inner_disk = gmsh.model.occ.addDisk(x0, y0, 0, r_inner, r_inner)
    vessel_wall, _ = gmsh.model.occ.cut([(2, outer_disk)], [(2, inner_disk)])

    # Coils geometry
    coil_rects = []
    for coil in coil_positions:
        x_min, x_max, y_min, y_max = coil
        rect = gmsh.model.occ.addRectangle(x_min, y_min, 0, x_max - x_min, y_max - y_min)
        coil_rects.append((2, rect))

    # Outer domain box
    domain_box = gmsh.model.occ.addRectangle(0, 0, 0, domain_w, domain_h)

    # Fragment all regions
    all_regions = [(2, domain_box), (2, inner_disk)] + vessel_wall + coil_rects
    gmsh.model.occ.fragment(all_regions, [])
    gmsh.model.occ.synchronize()

    # Classify physical regions
    surfaces = gmsh.model.occ.getEntities(dim=2)
    vacuum = []
    vessel = []
    coils = []
    air = []

    for tag in surfaces:
        com = gmsh.model.occ.getCenterOfMass(*tag)
        x, y = com[0], com[1]
        r = ((x - x0)**2 + (y - y0)**2)**0.5

        if r < r_inner - 1e-3:
            vacuum.append(tag[1])
        elif r_inner < r < r_outer + 1e-3:
            vessel.append(tag[1])
        elif any((x_min <= x <= x_max and y_min <= y <= y_max) for (x_min, x_max, y_min, y_max) in coil_positions):
            coils.append(tag[1])
        else:
            air.append(tag[1])

    gmsh.model.addPhysicalGroup(2, vacuum, tag=1)
    gmsh.model.setPhysicalName(2, 1, "Vacuum")
    gmsh.model.addPhysicalGroup(2, vessel, tag=2)
    gmsh.model.setPhysicalName(2, 2, "VesselWall")
    gmsh.model.addPhysicalGroup(2, coils, tag=3)
    gmsh.model.setPhysicalName(2, 3, "Coils")
    gmsh.model.addPhysicalGroup(2, air, tag=4)
    gmsh.model.setPhysicalName(2, 4, "Air")

    # Mesh options
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)

    # Finer mesh in vacuum region
    if vacuum:
        field_id = gmsh.model.mesh.field.add("Constant")
        gmsh.model.mesh.field.setNumbers(field_id, "SurfacesList", vacuum)
        gmsh.model.mesh.field.setNumber(field_id, "VIn", vacuum_mesh_size)
        gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

    gmsh.model.mesh.generate(2)
    gmsh.write(msh_path)
    gmsh.finalize()