import gmsh

def loop_from_points(points):
    """
    Create a closed loop from a list of points.
    pts: list of gmsh points ordered either clockwise or counterclockwise.
    Returns a tuple (lines, loop) where:
        - lines: list of line tags forming the loop
        - loop: counterclockwise oriented curve loop tag
    """

    lines = []
    for i in range(len(points) - 1):
        lines.append(gmsh.model.geo.addLine(points[i], points[i+1]))
    lines.append(gmsh.model.geo.addLine(points[-1], points[0]))

    loop = gmsh.model.geo.addCurveLoop(lines)

    return lines, loop


def check_counterclockwise(points):

    """
    Check if the loop is counterclockwise.
    points: list of points [(x1, y1), (x2, y2), ...]
    Returns True if the loop is counterclockwise, False otherwise.
    """

    area = 0.0 # Polygon area (times 2 since only the sign is needed)
    n = len(points)

    for i in range(n):
        j = (i + 1) % n
        x1, y1 = points[i]
        x2, y2 = points[j]
        area += x1 * y2 - x2 * y1

    # Positive area means counterclockwise orientation:
    return area > 0


def change_orientation(lines):
    """
    Defines a closed loop with reversed orientation.

    lines: list of line tags
    Returns a new curve loop tag with reversed orientation.
    """
    # Define 
    lines_rev = []
    for i in range(len(lines)):
        lines_rev.append(-lines[i])  # Invert the line orientation
    lines_rev = lines_rev[::-1] # Reverse the order of lines

    loop_rev = gmsh.model.geo.addCurveLoop(lines_rev)

    return loop_rev


def generate_mesh(params: dict, msh_path: str):
    """
    Use Gmsh Python API to build and save a Firedrake-compatible .msh mesh.
    params: dict with keys:
        - boundary: list of tuples defining the artificial boundary as a spline,
        - outer_wall: list of tuples defining the outer wall of the vacuum vessel,
        - inner_wall: list of tuples defining the inner wall of the vacuum vessel,
        - coils: list of tuples defining the coil geometries,
        - limiter_pts: list of tuples defining the limiter as a discrete set of points,
        - limiter_line: list of tuples defining the limiter as a continuos line (as a spline),
        - mesh_size_min: minimum mesh size (default 0.01),
        - mesh_size_max: maximum mesh size (default 0.05),
        - limiter_mesh_size: mesh size near the limiter (default 0.02),
        - dist_max: distance from the limiter up to which the limiter_mesh_size is applied (default 0.1).
    msh_path: output .msh file path
    """
    gmsh.initialize()
    gmsh.model.add("custom_tokamak")


    #--------------------------------------------------#
    #             DEFINE POINTS AND LINES              #
    #--------------------------------------------------#
    
    # Define artificial boundary closed loop:
    boundary = params.get("boundary", [(0.0,0.0),(1.0, 0.0),(1.0,1.0),(0.0,1.0)])
    if not check_counterclockwise(boundary):    # make sure points are ordered counterclockwise
        boundary = boundary[::-1]
    bdry_pts = [gmsh.model.geo.addPoint(float(xi), float(yi), 0) for xi, yi in boundary]

    bdry_lines, bdry_loop = loop_from_points(bdry_pts)

    # Define a close loop for each coil:
    coils = params["coils"]

    coils_lines = []
    coils_loops = []

    for coil in coils:
        if not check_counterclockwise(coil): # make sure points are ordered counterclockwise
            coil = coil[::-1]
        pts = [gmsh.model.geo.addPoint(float(xi), float(yi), 0) for xi, yi in coil]
        lines, loop = loop_from_points(pts)
        coils_lines.append(lines)
        coils_loops.append(loop)

    # Define a close loop for the walls of the vacuum vessel:
    outer_wall = params.get("outer_wall")
    inner_wall = params.get("inner_wall")

    # Make sure points are ordered counterclockwise:
    if not check_counterclockwise(outer_wall):
        outer_wall = outer_wall[::-1]
    if not check_counterclockwise(inner_wall):
        inner_wall = inner_wall[::-1]

    out_wll_pts = [gmsh.model.geo.addPoint(float(xi), float(yi), 0) for xi, yi in outer_wall]
    inn_wll_pts = [gmsh.model.geo.addPoint(float(xi), float(yi), 0) for xi, yi in inner_wall]

    out_wll_lines, out_wll_loop = loop_from_points(out_wll_pts)
    inn_wll_lines, inn_wll_loop = loop_from_points(inn_wll_pts)


    #--------------------------------------------------#
    #         DEFINE DOMAIN REGIONS AND LIMITER        #
    #--------------------------------------------------#

    # Define domain regions:
    coils_holes = []    # Invert the coils loops orientation to clockwise
    for loop in coils_loops:
        coils_holes.append(-loop)  # Invert the loop orientation to clockwise

    # Create the air region, vessel wall, and vacuum regions:
    air = gmsh.model.geo.addPlaneSurface([bdry_loop] + coils_holes + [-out_wll_loop])
    vessel_wall = gmsh.model.geo.addPlaneSurface([out_wll_loop, -inn_wll_loop])
    vacuum = gmsh.model.geo.addPlaneSurface([inn_wll_loop])

    # If limiter is a line, define it:
    if( params.get("limiter_line", None) is not None ):
        limiter = params["limiter_line"]
        if not check_counterclockwise(limiter): # make sure points are ordered counterclockwise
            limiter = limiter[::-1]
        lim_pts = [gmsh.model.geo.addPoint(float(xi), float(yi), 0) for xi, yi in limiter]
        lim_lines, lim_loop = loop_from_points(lim_pts)
        first_coil_tag = 5

    elif( params.get("limiter_pts", None) is not None ):
        lim_pts = params["limiter_pts"]
        n_pts = len(lim_pts)
        first_coil_tag = 4

        pt_list = []
        for i in range(n_pts):
            xi, yi = lim_pts[i]
            pt = [gmsh.model.geo.addPoint(float(xi), float(yi), 0)]
            pt_list.append(pt)

    # Define a region and a physical group for each coil:
    coils_list = []
    for i in range(len(coils_loops)):
        coil = gmsh.model.geo.addPlaneSurface([coils_loops[i]])
        coils_list.append(coil)

    # Synchronize
    gmsh.model.geo.synchronize()


    #--------------------------------------------------#
    #                 DEFINE PHYSICALS                 #
    #--------------------------------------------------#

    # Define physical gsurfaces for air, vessel wall, vacuum:
    gmsh.model.addPhysicalGroup(2, [air], tag=1)
    gmsh.model.setPhysicalName(2, 1, "Air Region")
    gmsh.model.addPhysicalGroup(2, [vessel_wall], tag=2)
    gmsh.model.setPhysicalName(2, 2, "Vessel Wall")
    gmsh.model.addPhysicalGroup(2, [vacuum], tag=3)
    gmsh.model.setPhysicalName(2, 3, "Vacuum Region")
    
    # Define physical line/points for the limiter:
    if( params.get("limiter_line", None) is not None ):
        # Embed into vacuum:
        gmsh.model.mesh.embed(1, [lim_lines], 2, vacuum)

        # Define physical line:
        gmsh.model.addPhysicalGroup(1, [lim_lines], tag=4)
        gmsh.model.setPhysicalName(1, 4, "Limiter")

    elif( params.get("limiter_pts", None) is not None ):
        i = 0
        for pt in pt_list:
            # Embed into vacuum:
            gmsh.model.mesh.embed(0, pt, 2, vacuum)

            # Define physical point:
            pt_tag = gmsh.model.addPhysicalGroup(0, pt)
            gmsh.model.setPhysicalName(0, pt_tag, "LimiterPoint_{}".format(i+1))

    # Define physical surfaces for the coils:
    i = 0
    for coil in coils_list:
        gmsh.model.addPhysicalGroup(2, [coil], tag = first_coil_tag + i)
        gmsh.model.setPhysicalName(2, first_coil_tag + i, "Coil_{}".format(i+1))
        i += 1

    #--------------------------------------------------#
    #                 MESH REFINEMENT                  #
    #--------------------------------------------------#

    mesh_size_min = params.get("mesh_size_min", 0.01)
    mesh_size_max = params.get("mesh_size_max", 0.05)
    limiter_mesh_size = params.get("limiter_mesh_size", 0.02)
    limiter_dist_max = params.get("limiter_dist_max", 0.0)

    # If the limiter is a closed loop, increase the refinment inside the limiter,
    # otherwise, refine around the limiter points:

    if params.get("limiter_line", None) is not None:
        # Refinement field for region inside the limiter:
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "CurvesList", lim_lines)
        gmsh.model.mesh.field.setNumber(distance, "NumPointsPerCurve", 100)

        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", limiter_mesh_size)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", mesh_size_max)
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", mesh_size_min)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", limiter_dist_max)    

    elif params.get("limiter_pts", None) is not None:
        # Refinement field for points around the limiter:
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "NodesList", [pt[0] for pt in pt_list])
        gmsh.model.mesh.field.setNumber(distance, "NumPointsPerCurve", 100)

        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", limiter_mesh_size)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", mesh_size_max)
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", limiter_dist_max)

    # Refine the mesh inside the coils:
    coils_mesh_size = params.get("coils_mesh_size", mesh_size_min)
    coils_dist_max = params.get("coils_dist_max", 0.0)
    coil_threshold_fields = []
    for i, lines in enumerate(coils_lines):

        # Create a distance field for each coil (using its boundary lines)
        field_dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(field_dist, "CurvesList", lines)
        gmsh.model.mesh.field.setNumber(field_dist, "NumPointsPerCurve", 50)
    
        # Create a threshold field for each coil
        field_thresh = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(field_thresh, "InField", field_dist)
        gmsh.model.mesh.field.setNumber(field_thresh, "SizeMin", coils_mesh_size)  # Finer mesh inside coil
        gmsh.model.mesh.field.setNumber(field_thresh, "SizeMax", mesh_size_max)   # Coarser mesh outside
        gmsh.model.mesh.field.setNumber(field_thresh, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(field_thresh, "DistMax", coils_dist_max)

        coil_threshold_fields.append(field_thresh)

    # Set background mesh field to combine all thresholds:
    all_thresholds = [threshold] + coil_threshold_fields
    final_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(final_field, "FieldsList", all_thresholds)

    gmsh.model.mesh.field.setAsBackgroundMesh(final_field)

    #--------------------------------------------------#
    #                 GENERATE MESH                    #
    #--------------------------------------------------#

    # Mesh options for Firedrake
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    #gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 2)

    # Generate and export
    gmsh.model.mesh.generate(2)
    gmsh.write(msh_path)
    gmsh.finalize()

    print("Mesh generated.")