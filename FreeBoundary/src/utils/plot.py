from firedrake import *
from firedrake.pyplot import triplot, tripcolor, tricontour
from src.utils.functions.mesh_tags import get_tags_for_plot
import matplotlib.pyplot as plt
import numpy as np


class Plot:
    '''
        Plots the mesh or the flux function provided together with
    '''

    def __init__(self, m, geometry):
        '''
            @brief Constructor.

            @param m: mesh
            @param geometry: tokamak geometry

            The constructor extract the coordinates on sets of nodes that lie on the tokamak
            structures (limiter, vessel walls, coils) to plot them with the flux.
        '''

        self.m = m
        V = FunctionSpace(m, "CG", 1)
        tags = get_tags_for_plot(geometry)

        # Get the coordinates of mesh nodes:
        coord_func = Function(VectorFunctionSpace(m, "CG", 1)).interpolate(as_vector(SpatialCoordinate(m)))
        coords = coord_func.dat.data_ro[:]
        
        # Get limiter nodes:
        limiter_idx = DirichletBC(V, 0.0, tags['limiter']).nodes
        lim = []
        for idx in limiter_idx:
            lim.append(coords[idx])
        self.limiter = np.array(lim)
    
        if geometry == 'custom':
            # Get inner wall points:
            inner_wall_idx = DirichletBC(V, 0.0, tags['inner_wall']).nodes
            iw = []
            for idx in inner_wall_idx:
                iw.append(coords[idx])
            self.inner_wall = np.array(iw)
            # Get outer wall points:
            outer_wall_idx = DirichletBC(V, 0.0, tags['outer_wall']).nodes
            ow = []
            for idx in outer_wall_idx:
                ow.append(coords[idx])
            self.outer_wall = np.array(ow)

            # ... coils missing ...

        elif geometry == "ITER":
            self.ITER_structures()

        else:
            raise ValueError(f"Geometry {geometry} not available.")


    def ITER_structures(self):
        '''
        Define the ITER wall structure. In ITER mesh the vessel is not embedded in the mesh,
        so it is not possible to recover its shape from the mesh file.
        Instead, the wall structure is defined manually by providing the coordinates of 
        some points of both the inner and external walls.
        '''

        # Points to plot the limiter
        lim = [
            [4.170213, -2.523364],[4.170213, 3.794393],[4.430851, 4.467290],[5.101064, 4.878505],
            [5.920213, 4.691589],[7.595745, 3.308411],[8.154255, 2.523364],[8.526596, 1.775701],
            [8.750000, 0.691589],[8.563830, -0.429907],[8.154255, -1.364486],[7.446809, -2.373832],
            [6.553191, -3.084112],[6.553191, -3.233645],[6.702128, -3.345794],[6.367021, -3.308411],
            [6.069149, -3.383178],[5.845745, -3.570093],[5.734043, -3.869159],[5.734043, -4.766355],
            [5.473404, -4.093458],[5.250000, -3.831776],[5.101064, -3.757009],[4.914894, -3.719626],
            [4.579787, -3.794393],[4.170213, -3.981308],[4.542553, -3.420561],[4.617021, -3.158879],
            [4.542553, -2.897196],[4.430851, -2.747664],[4.244681, -2.635514],[4.058511, -2.635514],
        ]

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

        # Coils:
        solenoid = np.array([[1.377660, 6.392523],[2.122340, 6.392523],[2.122340, -6.392523],[1.377660, -6.392523]])
        coil7 = np.array([[3.574468, 7.308411],[4.542553, 7.308411],[4.542553, 8.355140],[3.574468, 8.355140]])
        coil8 = np.array([[8.228723, 6.373832],[8.824468, 6.373832],[8.824468, 7.121495],[8.228723, 7.121495]])
        coil9 = np.array([[11.989362, 2.897196],[12.659574, 2.897196],[12.659574, 3.869159],[11.989362, 3.869159]])
        coil10 = np.array([[11.989362, -2.822430],[12.659574, -2.822430],[12.659574, -1.850467],[11.989362, -1.850467]])
        coil11 = np.array([[8.228723, -7.457944],[9.047872, -7.457944],[9.047872, -6.485981],[8.228723, -6.485981]])
        coil12 = np.array([[3.648936, -8.317757],[5.250000, -8.317757],[5.250000, -7.158879],[3.648936, -7.158879]])

        self.limiter = np.array(lim)
        self.inner_wall = np.array(iw)
        self.outer_wall = np.array(ow)
        self.coils_list = [solenoid, coil7, coil8, coil9, coil10, coil11, coil12]


    def mesh(self):
        '''
        @brief Plots the mesh
        '''
        print("Plotting mesh in file 'results/mesh_plot.png'...\n")
        fig, ax = plt.subplots()
        triplot(self.m, axes=ax)

        plt.title(r"Domain")
        plt.xlabel("r")
        plt.ylabel("z")
        plt.axis('equal')
        plt.savefig("./results/mesh_plot.png")
        plt.close()

    def flux(self,psi,psi0):
        '''
        @brief Plots the provided flux function psi.

        @param psi The flux function to plot.
        @param psi0 The flux contour line that identifies the plasdma boundary
        '''
        print("\nPlotting flux in file 'results/flux_plot.png'...\n")
        fig, ax = plt.subplots()

        # Plot flux and psi=psi0 contour line:
        fig.colorbar(tripcolor(psi,axes=ax))
        tricontour(psi, levels=[psi0], colors='red', linewidths=1.5, axes=ax)

        # Add limiter:
        ax.plot(self.limiter[:, 0], self.limiter[:, 1], 'b-', linewidth=1, label='Limiter')
        ax.plot([self.limiter[-1, 0], self.limiter[0, 0]], 
                [self.limiter[-1, 1], self.limiter[0, 1]], 'k-', linewidth=1)

        # Add tokamak walls:
        ax.plot(self.inner_wall[:, 0], self.inner_wall[:, 1], 'k-', linewidth=1, label='Inner wall')
        ax.plot(self.outer_wall[:, 0], self.outer_wall[:, 1], 'k-', linewidth=1, label='Outer wall')
        ax.plot([self.inner_wall[-1, 0], self.inner_wall[0, 0]], 
                [self.inner_wall[-1, 1], self.inner_wall[0, 1]], 'k-', linewidth=1)
        ax.plot([self.outer_wall[-1, 0], self.outer_wall[0, 0]], 
                [self.outer_wall[-1, 1], self.outer_wall[0, 1]], 'k-', linewidth=1)

        # Add coils:
        for coil in self.coils_list:
            ax.plot(coil[:, 0], coil[:, 1], 'k-', linewidth=1)
            ax.plot([coil[-1, 0], coil[0, 0]], 
                [coil[-1, 1], coil[0, 1]], 'k-', linewidth=1)

        plt.title(r"Flux")
        plt.xlabel("r")
        plt.ylabel("z")
        plt.axis('equal')
        plt.savefig("./results/flux_plot.png")
        plt.close()