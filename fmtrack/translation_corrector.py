import numpy as np
from pyearth import Earth
import matplotlib.pyplot as plt
import fmtrack
from tqdm import tqdm
from . import tracking
import fmtrack
import math

class TranslationCorrector:

    def __init__(self):
        pass

    # identifies 'safe' beads to use for translation correction by using
    # a bounding box (a rectangular prism) around the cell
    def bounding_box(self):

        # defines bounding box
        points_init = self.cell_init.points
        points_final = self.cell_final.points
        x_min = np.min([np.min(points_init[:,0]),np.min(points_final[:,0])]) - self.buffer_cell 
        x_max = np.max([np.max(points_init[:,0]),np.max(points_final[:,0])]) + self.buffer_cell
        y_min = np.min([np.min(points_init[:,1]),np.min(points_final[:,1])]) - self.buffer_cell
        y_max = np.max([np.max(points_init[:,1]),np.max(points_final[:,1])]) + self.buffer_cell
        z_min = np.min([np.min(points_init[:,2]),np.min(points_final[:,2])]) - self.buffer_cell
        z_max = np.max([np.max(points_init[:,2]),np.max(points_final[:,2])]) + self.buffer_cell

        # --> limit to points that aren't too close to the cell 
        X_safe = []; Y_safe = []; Z_safe = []; U_safe = []; V_safe = []; W_safe = [] 
        num_pts = len(self.U)
        for kk in range(0,num_pts):
            x_out = self.X[kk] < x_min or self.X[kk] > x_max
            y_out = self.Y[kk] < y_min or self.Y[kk] > y_max
            z_out = self.Z[kk] < z_min or self.Z[kk] > z_max
            if x_out or y_out or z_out:
                X_safe.append(self.X[kk])
                Y_safe.append(self.Y[kk])
                Z_safe.append(self.Z[kk])
                U_safe.append(self.U[kk])
                V_safe.append(self.V[kk])
                W_safe.append(self.W[kk])

        self.X_safe = np.asarray(X_safe); self.Y_safe = np.asarray(Y_safe); self.Z_safe = np.asarray(Z_safe)
        self.U_safe = np.asarray(U_safe); self.V_safe = np.asarray(V_safe); self.W_safe = np.asarray(W_safe)


    # identifies 'safe' beads to use for translation correction by using
    # only beads that are above a certain threshold from the cell boundary
    def distance_threshold(self):
        points = self.cell_init.points
        p1,p2,p3 = (points[:,0],points[:,1],points[:,2])

        b1,b2,b3 = (self.X,self.Y,self.Z)

        # Find shortest distance  between each initial bead and all the cell mesh vertices
        kk = 0
        new_beads = np.zeros((len(b1), 3))
        new_disps = np.zeros((len(b1), 3))
        temp = np.zeros(len(points))
        for bead in tqdm(range(0,len(b1)), desc='Computing Safe Beads', ascii=True):
            jj = 0
            for point in range(0,len(points)): 
                temp[jj] = math.sqrt((b1[bead]-p1[point])**2 + (b2[bead]-p2[point])**2 + (b3[bead]-p3[point])**2)
                jj = jj + 1
            
            min_dist = np.amin(temp)
            if min_dist > self.buffer_cell:
                new_beads[kk] = [b1[bead],b2[bead],b3[bead]]
                new_disps[kk] = [self.U[bead],self.V[bead],self.W[bead]]
                kk = kk + 1


        # Trim off the end zeros that were never filled
        new_beads = new_beads[:kk, :]
        new_disps = new_disps[:kk, :]
        self.X_safe, self.Y_safe, self.Z_safe = (new_beads[:,0],new_beads[:,1],new_beads[:,2])
        self.U_safe, self.V_safe, self.W_safe = (new_disps[:,0],new_disps[:,1],new_disps[:,2])

    # used only when there is no cell, or you want to include all of the beads in translation correction
    def no_cell(self):
        self.X_safe = self.X
        self.Y_safe = self.Y
        self.Z_safe = self.Z
        self.U_safe = self.U
        self.V_safe = self.V
        self.W_safe = self.W

    def create_model(self, beads_init, beads_final, cell_init, cell_final, buffer_cell, closest_no_conflict, use_box=True):

        # initialize objects
        self.beads_init = beads_init
        self.beads_final = beads_final
        self.cell_init = cell_init
        self.cell_final = cell_final
        self.buffer_cell = buffer_cell
        self.beads_init_new, self.beads_final_new = tracking.get_corrected_beads(beads_init, beads_final, closest_no_conflict)

        # get matched positions and displacements
        x_init, y_init, z_init = self.beads_init_new.get_xyz()
        self.X, self.Y, self.Z = self.beads_final_new.get_xyz()
        self.U = self.X - x_init
        self.V = self.Y - y_init
        self.W = self.Z - z_init

        # calculate safe beads
        if cell_init is None:
            self.no_cell()
        else:
            if use_box:
                self.bounding_box()
            else:
                self.distance_threshold()
        
        # --> fit MARS models 
        self.model_U = Earth(max_degree=2,max_terms=10)
        self.model_U.fit(self.Z_safe,self.U_safe)
        self.model_V = Earth(max_degree=2,max_terms=10)
        self.model_V.fit(self.Z_safe,self.V_safe)
        self.model_W = Earth(max_degree=2,max_terms=10)
        self.model_W.fit(self.Z_safe,self.W_safe)

    def correct_beads(self, beads):
        """Corrects beads using translation correction model

        Parameters
        ----------
        beads : FMBeads
            Input beads in 'final' state

        Returns
        ----------
        beads_new : FMBeads
            Corrected version of the input beads

        """

        x_pos, y_pos, z_pos = beads.get_xyz()

        # --> re-define Z 
        pred_U = self.model_U.predict(z_pos)
        pred_V = self.model_V.predict(z_pos)
        pred_W = self.model_W.predict(z_pos)
        
        # --> correct new bead positions 
        x_pos_new = x_pos.copy()
        y_pos_new = y_pos.copy()
        z_pos_new = z_pos.copy()
        for kk in range(0,len(x_pos)):
            x_pos_new[kk] = x_pos[kk] - pred_U[kk]
            y_pos_new[kk] = y_pos[kk] - pred_V[kk]
            z_pos_new[kk] = z_pos[kk] - pred_W[kk]

        beads_new = fmtrack.FMBeads()
        beads_new.load_from_positions(x_pos_new, y_pos_new, z_pos_new)

        return beads_new

    def correct_mesh(self, mesh):
        """Corrects mesh using translation correction model

        Parameters
        ----------
        mesh : FMMesh
            Input mesh in 'final' state

        Returns
        ----------
        mesh_new : FMMesh
            Corrected version of the input mesh

        """

        points = mesh.points

        # --> correct new cell position 
        pred_cell_0 = self.model_U.predict(points[:,2])
        pred_cell_1 = self.model_V.predict(points[:,2])
        pred_cell_2 = self.model_W.predict(points[:,2])
        
        points_new = np.zeros(points.shape)
        points_new[:,0] = points[:,0] - pred_cell_0
        points_new[:,1] = points[:,1] - pred_cell_1
        points_new[:,2] = points[:,2] - pred_cell_2

        mesh_new = fmtrack.FMMesh(points=points_new, faces=mesh.faces)
        
        return mesh_new
    
    def remove_spurious(self, beads_init_new, beads_final_new):
        """Removes spurious matches with displacements deviating greatly
        from translation correction predictions. Only removes 'far-field'
        beads, ones that are outside of the correction boundary (either a
        bounding box or distance threshold)

        Parameters
        ----------
        beads_init_new : FMBeads
            Matched beads in 
        beads_final_new : FMBeads

        Returns
        ----------
        mesh_new : FMMesh
            Corrected version of the input mesh

        """

    def create_figure(self):
        # --> plot MARS models 
        Z_line = np.linspace(np.min(self.Z),np.max(self.Z),100)
        pred_line_U = self.model_U.predict(Z_line)
        pred_line_V = self.model_V.predict(Z_line)
        pred_line_W = self.model_W.predict(Z_line)

        fig = plt.figure(figsize=(15,5))
        fig.tight_layout()
        ax1 = fig.add_subplot(1,3,1)
        ax1.plot(self.Z,self.U,'b.',label='x raw')
        ax1.plot(Z_line,pred_line_U,'k--',label='fit')
        ax1.set_xlabel('z position'); ax1.set_ylabel('displacement')
        ax1.legend(); ax1.set_title('x displacements')
        ax2 = fig.add_subplot(1,3,2)
        ax2.plot(self.Z,self.V,'r.',label='y raw')
        ax2.plot(Z_line,pred_line_V,'k--',label='fit')
        ax2.set_xlabel('z position'); ax2.set_ylabel('displacement')
        ax2.legend(); ax2.set_title('y displacements')
        ax3 = fig.add_subplot(1,3,3)
        ax3.plot(self.Z,self.W,'g.',label='z raw')
        ax3.plot(Z_line,pred_line_W,'k--',label='fit')
        ax3.set_xlabel('z position'); ax3.set_ylabel('displacement')
        ax3.legend(); ax3.set_title('z displacements')

        return fig
