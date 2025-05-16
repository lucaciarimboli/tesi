import gmsh
import meshio

gmsh.initialize()
gmsh.model.add("tokamak")

# Parameters
r_inner = 0.25
thickness = 0.03
x0, y0 = 0.7, 0.5
r_outer = r_inner + thickness

# === Step 1: Define Vessel (as two disks) ===
outer_disk = gmsh.model.occ.addDisk(x0, y0, 0, r_outer, r_outer)
inner_disk = gmsh.model.occ.addDisk(x0, y0, 0, r_inner, r_inner)

# === Step 2: Vessel wall = outer - inner ===
vessel_wall, _ = gmsh.model.occ.cut([(2, outer_disk)], [(2, inner_disk)])

# === Step 3: Define Coils ===
coil_rects = []
for coil in [(0.15, 0.25, 0.45, 0.55), (0.25, 0.35, 0.7, 0.8), (0.25, 0.35, 0.2, 0.3)]:
    x_min, x_max, y_min, y_max = coil
    rect = gmsh.model.occ.addRectangle(x_min, y_min, 0, x_max - x_min, y_max - y_min)
    coil_rects.append((2, rect))

# === Step 4: Define outer box (domain) ===
domain_box = gmsh.model.occ.addRectangle(0, 0, 0, 1.0, 1.0)

# === Step 5: Fragment everything together ===
all_regions = [(2, domain_box), (2, inner_disk)] + vessel_wall + coil_rects
gmsh.model.occ.fragment(all_regions, [])

# === Step 6: Synchronize geometry ===
gmsh.model.occ.synchronize()

# === Step 7: Classify physical regions ===

# Helper: get all surfaces
surfaces = gmsh.model.occ.getEntities(dim=2)

# Prepare containers
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
    elif any((x_min <= x <= x_max and y_min <= y <= y_max) for (x_min, x_max, y_min, y_max) in [
        (0.15, 0.25, 0.45, 0.55), (0.25, 0.35, 0.7, 0.8), (0.25, 0.35, 0.2, 0.3)
    ]):
        coils.append(tag[1])
    else:
        air.append(tag[1])

# Add physical groups
gmsh.model.addPhysicalGroup(2, vacuum, tag=1)
gmsh.model.setPhysicalName(2, 1, "Vacuum")

gmsh.model.addPhysicalGroup(2, vessel, tag=2)
gmsh.model.setPhysicalName(2, 2, "VesselWall")

gmsh.model.addPhysicalGroup(2, coils, tag=3)
gmsh.model.setPhysicalName(2, 3, "Coils")

gmsh.model.addPhysicalGroup(2, air, tag=4)
gmsh.model.setPhysicalName(2, 4, "Air")

# === Step 8: Mesh options ===
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.01)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.05)

# --- Add a mesh size field for the vacuum region (region 1) ---
# Get all surface tags for the vacuum region
vacuum_surfaces = vacuum  # already a list of surface tags

# Create a size field that is small on the vacuum region
field_id = gmsh.model.mesh.field.add("Constant")
gmsh.model.mesh.field.setNumbers(field_id, "SurfacesList", vacuum_surfaces)
gmsh.model.mesh.field.setNumber(field_id, "VIn", 0.02)  # Finer mesh in vacuum

# Set the background mesh field
gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

# --- End mesh size field ---

# Generate and save mesh
gmsh.model.mesh.generate(2)
gmsh.write("../Meshes/tokamak/tokamak_mesh.msh")
gmsh.finalize()