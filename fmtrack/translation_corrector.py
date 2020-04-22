import numpy as np
from pyearth import Earth
import matplotlib.pyplot as plt

class TranslationCorrector:

    def __init__(self):
        pass

    def create_model(self, beads_init, beads_final, cell_init, cell_final, buffer_cell, closest_no_conflict):
        if cell_init is not None:
            points_init = cell_init.points
            points_final = cell_final.points
        
            x_min = np.min([np.min(points_init[:,0]),np.min(points_final[:,0])]) - buffer_cell 
            x_max = np.max([np.max(points_init[:,0]),np.max(points_final[:,0])]) + buffer_cell
            y_min = np.min([np.min(points_init[:,1]),np.min(points_final[:,1])]) - buffer_cell
            y_max = np.max([np.max(points_init[:,1]),np.max(points_final[:,1])]) + buffer_cell
            z_min = np.min([np.min(points_init[:,2]),np.min(points_final[:,2])]) - buffer_cell
            z_max = np.max([np.max(points_init[:,2]),np.max(points_final[:,2])]) + buffer_cell

        x_pos, y_pos, z_pos = beads_init.get_xyz()
        x_pos_new, y_pos_new, z_pos_new = beads_init.get_xyz()
        
        num_pts = len(x_pos)
        X = []; Y = []; Z = []; U = []; V = []; W = [] 
        for kk in range(0,num_pts):
            idx = closest_no_conflict[kk]
            if idx < len(closest_no_conflict):
                U.append(x_pos_new[idx] - x_pos[kk])
                V.append(y_pos_new[idx] - y_pos[kk])
                W.append(z_pos_new[idx] - z_pos[kk])
                X.append(x_pos_new[idx]); Y.append(y_pos_new[idx]); Z.append(z_pos_new[idx])
        
        if cell_init is not None:
            # --> limit to points that aren't too close to the cell 
            X_safe = []; Y_safe = []; Z_safe = []; U_safe = []; V_safe = []; W_safe = [] 
            num_pts = len(U)
            for kk in range(0,num_pts):
                x_out = X[kk] < x_min or X[kk] > x_max
                y_out = Y[kk] < y_min or Y[kk] > y_max
                z_out = Z[kk] < z_min or Z[kk] > z_max
                if x_out or y_out or z_out:
                    X_safe.append(X[kk])
                    Y_safe.append(Y[kk])
                    Z_safe.append(Z[kk])
                    U_safe.append(U[kk])
                    V_safe.append(V[kk])
                    W_safe.append(W[kk])

            X_safe = np.asarray(X_safe); Y_safe = np.asarray(Y_safe); Z_safe = np.asarray(Z_safe)
            U_safe = np.asarray(U_safe); V_safe = np.asarray(V_safe); W_safe = np.asarray(W_safe)
        else:
            X_safe = X; Y_safe = Y; Z_safe = Z
            U_safe = U; V_safe = V; W_safe = W

        # saved for plotting
        self.U = U
        self.V = V
        self.W = W
        self.Z = Z
        
        # --> fit MARS models 
        self.model_U = Earth(max_degree=2,max_terms=10)
        self.model_U.fit(Z_safe,U_safe)
        self.model_V = Earth(max_degree=2,max_terms=10)
        self.model_V.fit(Z_safe,V_safe)
        self.model_W = Earth(max_degree=2,max_terms=10)
        self.model_W.fit(Z_safe,W_safe)

    def correct_beads(self, beads):
        x_pos, y_pos, z_pos = beads.get_xyz()

        # --> re-define Z 
        pred_U = self.model_U.predict(z_pos)
        pred_V = self.model_V.predict(z_pos)
        pred_W = self.model_W.predict(z_pos)
        
        # --> correct new bead positions 
        for kk in range(0,len(x_pos)):
            x_pos[kk] = x_pos[kk] - pred_U[kk] 
            y_pos[kk] = y_pos[kk] - pred_V[kk]
            z_pos[kk] = z_pos[kk] - pred_W[kk]

        beads.load_from_positions(x_pos, y_pos, z_pos)

        return beads

    def correct_mesh(self, mesh):
        points = mesh.points

        # --> correct new cell position 
        pred_cell_0 = self.model_U.predict(points[:,0])
        pred_cell_1 = self.model_V.predict(points[:,1])
        pred_cell_2 = self.model_W.predict(points[:,2])
        
        points_new = np.zeros(points.shape)
        points_new[:,0] = points[:,0] - pred_cell_0
        points_new[:,1] = points[:,1] - pred_cell_1
        points_new[:,2] = points[:,2] - pred_cell_2

        mesh_new = mesh
        mesh_new.points = points_new
        
        return mesh_new

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