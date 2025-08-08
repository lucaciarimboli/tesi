from firedrake import *
from firedrake.pyplot import triplot, tricontour, tricontourf
from src.utils.functions.mesh_tags import get_tags_for_plot
import matplotlib.pyplot as plt
import numpy as np


class Plot:
    '''
        Stores mesh and points, plots        
    '''

    def __init__(self, m, geometry, params = []):
        
        tags = get_tags_for_plot(geometry)

        self.limiter = DirichletBC(m, 0.0, tags['limiter']).nodes
        self.solenoid = DirichletBC(m, 0.0, tags['solenoid']).nodes
        self.coil7 = DirichletBC(m, 0.0, tags['coils'][0]).nodes
        self.coil8 = DirichletBC(m, 0.0, tags['coils'][1]).nodes
        self.coil9 = DirichletBC(m, 0.0, tags['coils'][2]).nodes
        self.coil10 = DirichletBC(m, 0.0, tags['coils'][3]).nodes
        self.coil11 = DirichletBC(m, 0.0, tags['coils'][4]).nodes
        self.coil12 = DirichletBC(m, 0.0, tags['coils'][5]).nodes

        if geometry == "custom":
            self.inner_wall = DirichletBC(m, 0.0, tags['inner_wall']).nodes
            self.outer_wall = DirichletBC(m, 0.0, tags['outer_wall']).nodes

        elif geometry == "ITER":
            self.set_ITER_wall_structure()


    def set_ITER_wall_structure(self):   

        # Points to plot inner wall
        iw = [
            [3.611702, 3.682243],[3.723404, 4.168224],[3.872340, 4.504673],[4.095745, 4.803738],
            [4.356383, 5.028037],[4.654255, 5.214953],[5.026596, 5.327103],[5.361702, 5.364486],
            [5.808511, 5.289720],[6.143617, 5.140187],[6.553191, 4.878505],[6.962766, 4.579439],
            [7.335106, 4.280374],[7.819149, 3.831776],[8.191489, 3.383178],[8.787234, 2.373832],
            [9.047872, 1.700935],[9.196809, 0.915888],[9.196809, 0.280374],[9.085106, -0.392523],
            [8.861702, -1.065421],[8.638298, -1.626168],[8.340426, -2.336449],[8.042553, -3.009346],
            [7.707447, -3.794393],[7.260638, -4.467290],[6.851064, -4.841121],[6.367021, -5.102804],
            [5.845745, -5.289720],[5.324468, -5.364486],[4.728723, -5.252336],[4.207447, -4.915888],
            [3.835106, -4.467290],[3.611702, -3.682243]
        ]
        # Points to plot outer wall
        ow = [
            [3.351064, 3.981308],[3.537234, 4.616822],[3.872340, 5.102804],[4.207447, 5.401869],
            [4.728723, 5.700935],[5.473404, 5.850467],[6.106383, 5.775701],[6.813830, 5.476636],
            [7.446809, 5.102804],[8.117021, 4.542056],[8.750000, 3.869159],[9.159574, 3.271028],
            [9.457447, 2.672897],[9.755319, 1.813084],[9.904255, 1.028037],[9.904255, 0.168224],
            [9.792553, -0.542056],[9.606383, -1.140187],[9.308511, -1.887850],[8.936170, -2.747664],
            [8.526596, -3.682243],[8.005319, -4.579439],[7.521277, -5.102804],[7.037234, -5.439252],
            [6.478723, -5.700935],[5.957447, -5.813084],[5.436170, -5.887850],[5.026596, -5.850467],
            [4.505319, -5.700935],[3.909574, -5.289720],[3.574468, -4.841121],[3.351064, -4.168224]
        ]
        self.inner_wall = np.array(iw)
        self.outer_wall = np.array(ow)


    def mesh(self):
        print("Plotting mesh in file 'results/mesh_plot.png'...")
        fig, ax = plt.subplots()
        triplot(self.mesh, axes=ax)
        #plt.scatter(*zip(*self.limiter), color='blue', label='Limiter Points')

        #if( self.params.get("geometry", "build") == "build" ):
            # Plot coils
        #    for coil in self.params["coils_adapted_format"]:
        #        x_min, x_max, y_min, y_max = coil
        #        plt.plot([x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], color='red', label='Coil Edges')

        plt.title(r"Domain")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis('equal')
        plt.savefig("./results/mesh_plot.png")
        plt.close()