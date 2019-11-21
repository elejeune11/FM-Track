import numpy as np 
import os
from . import tracking
from . import fmplot
from . import post_process
from . import fmbeads
from . import fmmesh
import pickle

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)
from sklearn.neighbors import KernelDensity
from sklearn import preprocessing

class FMTracker:

	def __init__(self,cell_init=None,cell_final=None,beads_init=None,beads_final=None):
        
		self.cell_init = cell_init
		self.cell_final = cell_final

		self.beads_init = beads_init
		self.beads_final = beads_final

		self.beads_init_new = None
		self.beads_final_new = None

		self.closest_no_conflict = None

		self.print_progress = True

		# For get_tracking_params()
		self.num_feat = 5
		self.num_nearest = 15
		self.buffer_cell = 0
		self.track_type = 2 # type 1 will NOT perform translation correction, type 2 will

	def save(self,filename):
		pickle.dump(self, open(filename,'wb'))

	def run_tracking(self):
		num_feat, num_nearest, buffer_cell, track_type  = self.get_tracking_params()
		self.closest_no_conflict, self.idx_ignored, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, cell_final_new, self.mars_figure = \
			tracking.track_main_call(track_type,self.beads_init.points, self.beads_final.points, self.cell_init, self.cell_final, num_feat, num_nearest, buffer_cell, self.print_progress)
		self.cell_final_new = cell_final_new
		self.beads_init_new = fmbeads.FMBeads(np.transpose(np.vstack((x_pos,y_pos,z_pos))))
		self.beads_final_new = fmbeads.FMBeads(np.transpose(np.vstack((x_pos_new,y_pos_new,z_pos_new))))

	def save_res(self,folder,label_uncorrected):
		x_pos = self.beads_init_new.points[:,0]
		y_pos = self.beads_init_new.points[:,1]
		z_pos = self.beads_init_new.points[:,2]

		x_pos_new = self.beads_final_new.points[:,0]
		y_pos_new = self.beads_final_new.points[:,1]
		z_pos_new = self.beads_final_new.points[:,2]

		tracking.save_res(folder, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, self.closest_no_conflict, label_uncorrected)

	def load_res(self,bead_folder,cell_init_name,cell_final_name,label_uncorrected):
		x_pos, y_pos, z_pos, U, V, W = tracking.load_res(bead_folder,label_uncorrected)
		x_pos_new = x_pos + U
		y_pos_new = y_pos + V
		z_pos_new = z_pos + W

		self.beads_init_new = fmbeads.FMBeads()
		self.beads_final_new = fmbeads.FMBeads()

		self.beads_init_new.points = np.transpose(np.vstack((x_pos,y_pos,z_pos)))
		self.beads_final_new.points = np.transpose(np.vstack((x_pos_new,y_pos_new,z_pos_new)))

		self.cell_init = fmmesh.FMMesh()
		self.cell_init.import_native_files(cell_init_name)
		self.cell_final = fmmesh.FMMesh()
		self.cell_final.import_native_files(cell_final_name)

	def save_mars_figure(self,filename):
		self.mars_figure.savefig(filename)

	def plot(self):
		plotter = fmplot.FMPlot(self)
		plotter.plot()

	def create_gp_model(self):
		X = self.beads_init_new.points[:,0]
		Y = self.beads_init_new.points[:,1]
		Z = self.beads_init_new.points[:,2]
		U = self.beads_final_new.points[:,0] - Xf
		V = self.beads_final_new.points[:,1] - Y
		W = self.beads_final_new.points[:,2] - Z

		X, Y, Z, U, V, W = threshold(X, Y, Z, U, V, W, 3, False)
		self.gp_U, _ = create_gp_model(X,Y,Z,U)
		self.gp_V, _ = create_gp_model(X,Y,Z,V)
		self.gp_W, self.scaler = create_gp_model(X,Y,Z,W)

	def save_gp_model(self,foldername):
		pickle.dump(self.gp_U, open(os.path.join(foldername,'gp_U.sav'),'wb'))
		pickle.dump(self.gp_V, open(os.path.join(foldername,'gp_V.sav'),'wb'))
		pickle.dump(self.gp_W, open(os.path.join(foldername,'gp_W.sav'),'wb'))
		pickle.dump(self.scaler,open(os.path.join(foldername,'scaler.sav'),'wb'))

	def load_gp_model(self,foldername):
		pickle.load(self.gp_U, open(os.path.join(foldername,'gp_U.sav'),'rb'))
		pickle.load(self.gp_V, open(os.path.join(foldername,'gp_V.sav'),'rb'))
		pickle.load(self.gp_W, open(os.path.join(foldername,'gp_W.sav'),'rb'))
		pickle.load(self.scaler,open(os.path.join(foldername,'scaler.sav'),'rb'))

	def get_tracking_params(self):
		return self.num_feat, self.num_nearest, self.buffer_cell, self.track_type 

def load_fmtrack(filename):
    return pickle.load(open(filename, 'rb'))

def create_gp_model(X,Y,Z,QoI):
	num_pts = X.shape[0]
	X_train_unscale = np.zeros((num_pts,3))
	X_train_unscale[:,0] = X
	X_train_unscale[:,1] = Y
	X_train_unscale[:,2] = Z 
	scaler = preprocessing.StandardScaler().fit(X_train_unscale)
	X_train = scaler.transform(X_train_unscale)
	kernel = RationalQuadratic()
	gp = GaussianProcessRegressor(kernel=kernel)
	gp.fit(X_train, QoI)
	return gp , scaler

def threshold(X, Y, Z, U, V, W, thresh, above):
	XR = np.array([])
	YR = np.array([])
	ZR = np.array([])
	UR = np.array([])
	VR = np.array([])
	WR = np.array([])

	for i in range(X.shape[0]):
		if above:
			if np.linalg.norm(np.array([U[i],V[i],W[i]])) >= thresh:
				XR = np.append(XR, X[i])
				YR = np.append(YR, Y[i])
				ZR = np.append(ZR, Z[i])
				UR = np.append(UR, U[i])
				VR = np.append(VR, V[i])
				WR = np.append(WR, W[i])
		else:
			if np.linalg.norm(np.array([U[i],V[i],W[i]])) <= thresh:
				XR = np.append(XR, X[i])
				YR = np.append(YR, Y[i])
				ZR = np.append(ZR, Z[i])
				UR = np.append(UR, U[i])
				VR = np.append(VR, V[i])
				WR = np.append(WR, W[i])
	return XR, YR, ZR, UR, VR, WR