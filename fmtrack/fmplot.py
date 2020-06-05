import pickle
import pyvista
import numpy as np
import matplotlib.pyplot as plt
from . import post_process
from . import fmmesh
from . import fmtracker
import math
from pathlib import Path

class FMPlot:

    def __init__(self, *args, **kwargs):

        self.plot_at_bead_locations = True
        self.plot_displacement_vectors = False
        self.plot_principal_axes = False
        self.plot_principal_strains = False
        self.plot_jacobian = False

        # plot the cells
        self.plot_cell_init = False
        self.plot_cell_final = False

        # other parameters
        self.num_pts = 20
        self.eigen_num = 1    # 1 for largest eigenvalue, 3 for smallest
        self.color_by = 'cell_center'
        self.uniform_deformation_tensor = None
        self.nonuniform_deformation_tensor = None
        self.jacobians = None
        self.principal_axes = None
        self.cell_init = None
        self.cell_final = None
        self.cmap = plt.cm.get_cmap('coolwarm')
        self.text_color = 'white'
        self.scale_by_magnitude = False
        self.threshold_above = False
        self.threshold_below = False

        self.figtype_list = ['.png']

        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None
        self.z_min = None
        self.z_max = None
        self.X_DIM = None
        self.Y_DIM = None
        self.Z_DIM = None
        self.neigh_score = None
        self.mag_list = None
        self.dist_from_cell = None
        self.dist_from_edge = None
        self.dir_score = None

        self.num_feat = None
        if len(args) == 1:
            if isinstance(args[0], fmtracker.FMTracker):
                self._from_fmtracker(args[0])

    def _from_fmtracker(self,tracker):
        self.add_cell_init(tracker.cell_init)
        self.add_cell_final(tracker.cell_final)
        if tracker.beads_init_nonspurious is None:
            self.add_bead_positions(tracker.beads_init_new.points)
            self.add_bead_displacements(tracker.beads_final_new.points - tracker.beads_init_new.points)
        else:
            self.add_bead_positions(tracker.beads_init_nonspurious.points)
            self.add_bead_displacements(tracker.beads_final_nonspurious.points - tracker.beads_init_nonspurious.points)
            
        self.plot_displacement_vectors = True
        self.num_feat = tracker.num_feat
        if tracker.mars_model is not None:
            self.translation_correction = tracker.mars_model.create_figure()
        else:
            self.translation_correction = None

    def save(self,filename):
        self.plotter = None
        pickle.dump(self, open(filename,'wb'))

    def _before_any_2d_plot(self, recalculate=False):
        if self.num_feat is None:
            raise Exception('Must specify FMPlot.num_feat parameter before graphing')

        X = self.bead_positions[:,0]
        Y = self.bead_positions[:,1]
        Z = self.bead_positions[:,2]
        U = self.bead_displacements[:,0]
        V = self.bead_displacements[:,1]
        W = self.bead_displacements[:,2]
        if self.neigh_score is None or recalculate:
            self.neigh_score = post_process.color_point_neighbor_similarity(X, Y, Z, U, V, W, self.num_feat)
        if self.dir_score is None or self.dist_from_cell is None or self.mag_list is None or recalculate:
            self.dir_score, self.dist_from_cell, self.mag_list = post_process.color_point_direction(X, Y, Z, U, V, W, self.cell_init)
        if self.dist_from_edge is None or recalculate:
            self.dist_from_edge = post_process.compute_dist_from_edge(X, Y, Z, self.X_DIM, self.Y_DIM, self.Z_DIM)

    def save_plot_summary(self,filename):
        self._before_any_2d_plot()
        X = self.bead_positions[:,0]
        Y = self.bead_positions[:,1]
        Z = self.bead_positions[:,2]
        U = self.bead_displacements[:,0]
        V = self.bead_displacements[:,1]
        W = self.bead_displacements[:,2]
        post_process.plot_all(filename,self.dir_score,self.neigh_score,self.dist_from_edge,self.dist_from_cell,self.mag_list,\
			X,Y,Z,U,V,W,self.cell_init,self.cell_final,self.X_DIM,self.Y_DIM,self.Z_DIM,self.figtype_list)
    
    def save_plot_only_cells(self,filename):
        self._before_any_2d_plot()
        post_process.plot_only_cells(filename,self.cell_init,self.cell_final,self.X_DIM,self.Y_DIM,self.Z_DIM,self.figtype_list)

    def save_plot_only_scores(self,filename):
        self._before_any_2d_plot()
        post_process.plot_only_scores(filename,self.neigh_score,self.dir_score,self.figtype_list)

    def save_plot_only_slice(self,filename):
        self._before_any_2d_plot()
        X = self.bead_positions[:,0]
        Y = self.bead_positions[:,1]
        Z = self.bead_positions[:,2]
        U = self.bead_displacements[:,0]
        V = self.bead_displacements[:,1]
        W = self.bead_displacements[:,2]
        post_process.plot_only_slice(filename,self.dir_score,X,Y,Z,U,V,W,self.cell_init,self.cell_final,self.X_DIM,self.Y_DIM,self.Z_DIM,self.figtype_list)

    def save_plot_only_distance(self,filename):
        self._before_any_2d_plot()
        post_process.plot_only_distance(filename,self.cell_init,self.dist_from_edge,self.dist_from_cell,self.mag_list,self.figtype_list)

    def save_plot_only_translation_correction(self,filename):
        if self.translation_correction is not None:
            self.translation_correction.savefig(filename)

    def save_vtk_files(self,folder, cell_init_name='cell_init', cell_final_name='cell_final', arrows_name='arrows'):
        point_cloud = pyvista.PolyData(self.bead_positions)
        point_cloud["dot(cell normal, displacement)"] = self.dir_score
        point_cloud['vectors'] = self.bead_displacements
        geom = pyvista.Arrow()
        arrows = point_cloud.glyph(orient='vectors', scale=False, factor=5.0,geom=geom)

        mesh_init = pyvista.PolyData(self.cell_init.points)
        mesh_final = pyvista.PolyData(self.cell_final.points)

        mesh_init.save(Path(folder).joinpath(cell_init_name).with_suffix('.vtk'))
        mesh_final.save(Path(folder).joinpath(cell_final_name).with_suffix('.vtk'))
        arrows.save(Path(folder).joinpath(arrows_name).with_suffix('.vtk'))

    def save_native_plots(self,folder):
        path = Path(folder).joinpath('Post_proc_summary')
        self.save_plot_summary(path)
        path = Path(folder).joinpath('Cell_plots_3D')
        self.save_plot_only_cells(path)
        path = Path(folder).joinpath('Score_plots')
        self.save_plot_only_scores(path)
        path = Path(folder).joinpath('Bead_disp_slice')
        self.save_plot_only_slice(path)
        path = Path(folder).joinpath('Disp_wrt_dist')
        self.save_plot_only_distance(path)
        path = Path(folder).joinpath('Translation_correction')
        self.save_plot_only_translation_correction(path)

    def add_cell_init(self,cell):
        self.plot_cell_init = True
        self.cell_init = cell

    def add_cell_final(self,cell):
        self.plot_cell_final = True
        self.cell_final = cell

    def add_gpr_models(self,gp_U,gp_V,gp_W,scaler):
        self.gp_U = gp_U
        self.gp_V = gp_V
        self.gp_W = gp_W
        self.scaler = scaler

    def create_gp_model(self):
        X = self.bead_positions[:,0]
        Y = self.bead_positions[:,1]
        Z = self.bead_positions[:,2]
        U = self.bead_displacements[:,0]
        V = self.bead_displacements[:,1]
        W = self.bead_displacements[:,2]

        X, Y, Z, U, V, W = fmtracker.threshold(X, Y, Z, U, V, W, 3, False)
        self.gp_U, _ = fmtracker.create_gp_model(X,Y,Z,U)
        self.gp_V, _ = fmtracker.create_gp_model(X,Y,Z,V)
        self.gp_W, self.scaler = fmtracker.create_gp_model(X,Y,Z,W)

    def add_bead_positions(self,positions):
        self.bead_positions = positions
        self.x_min = np.amin(positions[:,0])
        self.x_max = np.amax(positions[:,0])
        self.y_min = np.amin(positions[:,1])
        self.y_max = np.amax(positions[:,1])
        self.z_min = np.amin(positions[:,2])
        self.z_max = np.amax(positions[:,2])

        if self.X_DIM is None:
            self.X_DIM = abs(self.x_max - self.x_min)
        if self.Y_DIM is None:
            self.Y_DIM = abs(self.y_max - self.y_min)
        if self.Z_DIM is None:
            self.Z_DIM = abs(self.z_max - self.z_min)

    def add_bead_displacements(self,displacements):
        self.bead_displacements = displacements

    def clear(self):
        self.remove_cell_init()
        self.remove_cell_final()
        self.remove_gpr_models()
        self.remove_bead_positions()
        self.remove_bead_displacements()

    def remove_cell_init(self):
        self.cell_init = None
        self.plot_cell_init = False

    def remove_cell_final(self):
        self.cell_final = None
        self.plot_cell_final = None

    def remove_gpr_models(self):
        self.gp_U = None
        self.gp_V = None
        self.gp_W = None
        self.scaler = None

    def remove_bead_positions(self):
        self.bead_positions = None
        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None
        self.z_min = None
        self.z_max = None

    def remove_bead_displacements(self):
        self.bead_displacements = None

    def set_min_max(self,x_min,x_max,y_min,y_max,z_min,z_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = x_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max

    def plot(self,title=None):
        self.plotter = pyvista.Plotter()
        if self.plot_cell_init:
            self._plot_cell_init()
        if self.plot_cell_final:
            self._plot_cell_final()
        if self.plot_displacement_vectors:
            self._plot_displacement_vectors()

        if self.plot_principal_axes:
            self._plot_principal_axes()
        if self.plot_jacobian:
            self._plot_jacobian()
        self.plotter.show_grid(color=self.text_color)
        self.plotter.view_xy(negative=True)
        if title is not None:
            self.plotter.show(title=title)
        else:
            self.plotter.show()

    def _plot_cell_init(self):
        self.plotter.add_mesh(self.cell_init.points,color='maroon')

    def _plot_cell_final(self):
        self.plotter.add_mesh(self.cell_final.points,color='maroon')

    def _plot_displacement_vectors(self):
        XYZ = self.bead_positions
        UVW = self.bead_displacements

        if self.color_by == 'cell_center' and self.threshold_above:
            X = np.array([])
            Y = np.array([])
            Z = np.array([])
            U = np.array([])
            V = np.array([])
            W = np.array([])
            dir_score = np.array([])
            score = get_color(XYZ, UVW, self.cell_init.center)
            for i in range(score.shape[0]):
                if score[i] > 0:

                    X = np.append(X, XYZ[i,0])
                    Y = np.append(Y, XYZ[i,1])
                    Z = np.append(Z, XYZ[i,2])
                    U = np.append(U, UVW[i,0])
                    V = np.append(V, UVW[i,1])
                    W = np.append(W, UVW[i,2])
                    dir_score = np.append(dir_score,score[i])

            XYZ = np.transpose(np.vstack((X,Y,Z)))
            UVW = np.transpose(np.vstack((U,V,W)))

        if self.color_by == 'cell_center' and self.threshold_below:
            X = np.array([])
            Y = np.array([])
            Z = np.array([])
            U = np.array([])
            V = np.array([])
            W = np.array([])
            dir_score = np.array([])
            score = get_color(XYZ, UVW, self.cell_init.center)
            for i in range(score.shape[0]):
                if score[i] < 0:

                    X = np.append(X, XYZ[i,0])
                    Y = np.append(Y, XYZ[i,1])
                    Z = np.append(Z, XYZ[i,2])
                    U = np.append(U, UVW[i,0])
                    V = np.append(V, UVW[i,1])
                    W = np.append(W, UVW[i,2])
                    dir_score = np.append(dir_score,score[i])

            XYZ = np.transpose(np.vstack((X,Y,Z)))
            UVW = np.transpose(np.vstack((U,V,W)))



            #dir_score = get_color(XYZ, UVW, self.cell_init.center)
            #X, Y, Z, _, _, _ = post_process.threshold(XYZ[:,0],XYZ[:,1],XYZ[:,2],dir_score[:,0],dir_score[:,1],dir_score[:,2],0,True)
            #XYZ = np.transpose(np.vstack((X,Y,Z)))
            #U, V, W, d1, d2, d3 = post_process.threshold(UVW[:,0],UVW[:,1],UVW[:,2],dir_score[:,0],dir_score[:,1],dir_score[:,2],0,True)
            #UVW = np.transpose(np.vstack((U,V,W)))
            #dir_score = np.transpose(np.vstack((d1,d2,d3)))




        point_cloud = pyvista.PolyData(XYZ)
        if self.cell_init is not None and self.color_by == 'dot':
            dir_score, _, _ = post_process.color_point_direction(XYZ[:,0], XYZ[:,1], XYZ[:,2], UVW[:,0], UVW[:,1], UVW[:,2], self.cell_init)
            point_cloud['Dot Product of Cell Normal and Displacement Vectors'] = dir_score
        elif self.color_by == 'cell_center':
            dir_score = get_color(XYZ, UVW, self.cell_init.center)
            point_cloud['dir_score'] = dir_score
        elif self.color_by == 'magnitude':
            point_cloud['dir_score'] = np.linalg.norm(UVW, axis=1)
        point_cloud['vectors'] = UVW
        geom = pyvista.Arrow()
        arrows = point_cloud.glyph(orient='vectors', scale=self.scale_by_magnitude, factor=5.0,geom=geom)

        self.plotter.add_mesh(arrows, cmap=self.cmap)
        self._create_scalar_bar()
            

    def _plot_jacobian(self):
        x, y, z = self._generate_uniform_grid()
        if self.uniform_deformation_tensor is None:
            self.uniform_deformation_tensor = post_process.calculate_deformation_tensor(x,y,z,self.gp_U,self.gp_V,self.gp_W,self.scaler)

        self.jacobians = np.array([])

        for F in self.uniform_deformation_tensor:
            self.jacobians = np.append(self.jacobians, abs(np.linalg.det(F)))

        grid = pyvista.UniformGrid()
        grid.dimensions = (self.num_pts,self.num_pts,self.num_pts)
        x_space = (self.x_max - self.x_min)/self.num_pts
        y_space = (self.y_max - self.y_min)/self.num_pts
        z_space = (self.z_max - self.z_min)/self.num_pts
        grid.spacing = (x_space, y_space, z_space) 
        grid.point_arrays["values"] = self.jacobians

        self.plotter.add_mesh_slice_orthogonal(grid, cmap=self.cmap)
        self.plotter.remove_scalar_bar()
        self.plotter.add_scalar_bar('Deformation Tensor Jacobian', title_font_size=20, label_font_size=15, color=self.text_color)

    def _plot_principal_axes(self):
        if self.plot_at_bead_locations:
            x = self.bead_positions[:,0]
            y = self.bead_positions[:,1]
            z = self.bead_positions[:,2]
            self.nonuniform_deformation_tensor = post_process.calculate_deformation_tensor(x,y,z,self.gp_U,self.gp_V,self.gp_W,self.scaler)
            deformation_tensor = self.nonuniform_deformation_tensor
        else:
            x, y, z = self._generate_uniform_grid()
            if self.uniform_deformation_tensor is None:
                self.uniform_deformation_tensor = post_process.calculate_deformation_tensor(x,y,z,self.gp_U,self.gp_V,self.gp_W,self.scaler)
            deformation_tensor = self.uniform_deformation_tensor

        X, Y, Z = np.meshgrid(x,y,z)
        X = X.flatten()
        Y = Y.flatten()
        Z = Z.flatten()
        XYZ = np.transpose(np.vstack((X,Y,Z)))

        self.principal_axes = np.empty((0,3), float)
        self.eigenvalues = np.array([])

        for F in deformation_tensor:
            w,v = np.linalg.eig(F)
            if self.eigen_num == 1:
                i = np.argmax(np.abs(w))
            elif self.eigen_num == 3:
                i = np.argmin(np.abs(w))
            self.principal_axes = np.vstack((self.principal_axes,np.transpose(v[:,i].real)))
            self.eigenvalues = np.append(self.eigenvalues, w[i].real)

        point_cloud = pyvista.PolyData(XYZ)
        point_cloud['vectors'] = self.principal_axes
        if self.cell_init is not None and self.color_by == 'dot':
            UVW = self.principal_axes
            dir_score, _, _ = post_process.color_point_direction(XYZ[:,0], XYZ[:,1], XYZ[:,2], UVW[:,0], UVW[:,1], UVW[:,2], self.cell_init)
            point_cloud['Dot Product of Cell Normal and Displacement Vectors'] = dir_score
        elif self.color_by == 'eigenvalue':
            point_cloud['Eigenvalue'] = self.eigenvalues
    
        geom = pyvista.Arrow()
        arrows = point_cloud.glyph(orient='vectors', scale=True, factor=1.0,geom=geom,tolerance=0.0)
        self.plotter.add_mesh(arrows, cmap=self.cmap)
        self._create_scalar_bar()

    def _generate_uniform_grid(self):
        if self.bead_positions is None:
            raise Exception('must set FMPlot.bead_positions or call FMPlot.set_min_max to set dimensions')

        x = np.linspace(self.x_min,self.x_max,self.num_pts)
        y = np.linspace(self.y_min,self.y_max,self.num_pts)
        z = np.linspace(self.z_min,self.z_max,self.num_pts)

        return x, y, z

    def _create_scalar_bar(self):
        self.plotter.remove_scalar_bar()
        if self.color_by == 'dot':
            self.plotter.add_scalar_bar('Dot(Cell Normal, Vector)', title_font_size=20, label_font_size=15, position_y=0.05, color=self.text_color)
        elif self.color_by == 'eigenvalue':
            self.plotter.add_scalar_bar('Eigenvalues of Principal Axes', title_font_size=20, label_font_size=15, position_y=0.05, color=self.text_color)
        elif self.color_by == 'cell_center':
            self.plotter.add_scalar_bar('Dot(Cell Center to Bead Unit Vector, Displacement)', title_font_size=20, label_font_size=15, position_y=0.02, position_x=0.8, color=self.text_color)
        elif self.color_by == 'magnitude':
            self.plotter.add_scalar_bar('Magnitude of Displacement', title_font_size=20, label_font_size=15, position_y=0.05, color=self.text_color)

    def save_screenshot(self,viewpoint):
        self.plotter.view_xy()
        self.plotter.screenshot('/Users/jakesansom/Desktop/image.pdf')

def load_fmplot(filename):
    return pickle.load(open(filename, 'rb'))

def get_color(XYZ, UVW, center):
    to_center = np.empty((0,3), float)
    for elem in XYZ:
        vec = (elem-center)/np.linalg.norm(elem-center)
        to_center = np.vstack((to_center,vec))

    dir_score = np.array([])
    for i in range(UVW.shape[0]):
        vec = UVW[i]#/np.linalg.norm(UVW[i])
        dir_score = np.append(dir_score, np.dot(vec,to_center[i,:]))

    dir_score = dir_score #/ np.nanmax(np.abs(dir_score))

    return dir_score
