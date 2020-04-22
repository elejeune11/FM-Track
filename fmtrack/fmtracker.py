import numpy as np 
import os
from . import tracking
from . import fmplot
from . import post_process
from . import fmbeads
from . import fmmesh
from . import translation_corrector
import pickle
from tqdm import tqdm

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

		self.num_feat = 5
		self.num_nearest = 15
		self.buffer_cell = 30
		self.track_type = 2 # type 1 will NOT perform translation correction, type 2 will

	def save(self,filename):
		# saves entire FMTracker object by pickling it
		pickle.dump(self, open(filename,'wb'))

	def save_res(self,folder,label_uncorrected):
		# saves bead positions in initial state along with displacements as text files
		tracking.save_res(folder, self.beads_init, self.beads_final, label_uncorrected)

	def load_res(self,bead_folder,cell_init_name,cell_final_name,label_uncorrected):
		# loads bead positions and cell meshes from native text file formats
		self.beads_init_new, self.beads_final_new = tracking.load_res(bead_folder,label_uncorrected)

		self.cell_init = fmmesh.FMMesh()
		self.cell_init.import_native_files(cell_init_name)
		self.cell_final = fmmesh.FMMesh()
		self.cell_final.import_native_files(cell_final_name)

	def save_mars_figure(self,filename):
		# saves the translation correction plot
		fig = self.mars_model.create_figure()
		fig.savefig(filename)

	def plot(self):
		# creates a 3D plot of the tracking results
		plotter = fmplot.FMPlot(self)
		plotter.plot()

	def create_gp_model(self):
		# creates a Gaussian process model
		X = self.beads_init_new.points[:,0]
		Y = self.beads_init_new.points[:,1]
		Z = self.beads_init_new.points[:,2]
		U = self.beads_final_new.points[:,0] - X
		V = self.beads_final_new.points[:,1] - Y
		W = self.beads_final_new.points[:,2] - Z

		X, Y, Z, U, V, W = threshold(X, Y, Z, U, V, W, 3, False)
		self.gp_U, _ = create_gp_model(X,Y,Z,U)
		self.gp_V, _ = create_gp_model(X,Y,Z,V)
		self.gp_W, self.scaler = create_gp_model(X,Y,Z,W)

	def save_gp_model(self,foldername):
		# saves the Gaussian process models to the specified folder
		pickle.dump(self.gp_U, open(os.path.join(foldername,'gp_U.sav'),'wb'))
		pickle.dump(self.gp_V, open(os.path.join(foldername,'gp_V.sav'),'wb'))
		pickle.dump(self.gp_W, open(os.path.join(foldername,'gp_W.sav'),'wb'))
		pickle.dump(self.scaler,open(os.path.join(foldername,'scaler.sav'),'wb'))

	def load_gp_model(self,foldername):
		# loads the Gaussian process models from the specified folder
		pickle.load(self.gp_U, open(os.path.join(foldername,'gp_U.sav'),'rb'))
		pickle.load(self.gp_V, open(os.path.join(foldername,'gp_V.sav'),'rb'))
		pickle.load(self.gp_W, open(os.path.join(foldername,'gp_W.sav'),'rb'))
		pickle.load(self.scaler,open(os.path.join(foldername,'scaler.sav'),'rb'))

	def run_tracking(self):
		"""Main tracking function. If self.track_type is 2, run_tracking() runs the more robust tracking algorithm,
		correcting for microscope translation. If self.track_type is 1, run_tracking() does not correct for microscope
		translation.

		"""
			
		if self.track_type == 1:
			closest_no_conflict = self.two_way_track()
		elif self.track_type == 2:
			closest_no_conflict = self.track_correct_track()

		self.beads_init_new, self.beads_final_new = tracking.get_corrected_beads(self.beads_init, self.beads_final, closest_no_conflict)

	def track_correct_track(self):
		"""Performs a round of two-way tracking, followed by translation correction, followed by a second round of two-way tracking

		Returns
		----------
		closest_no_conflict : np.array
			One dimensional NumPy array. Each element closest_no_conflict[i] contains the index of the bead in beads_final
			correlating to bead i in beads_init. If no corresponding bead was found, closest_no_conflict[i] = 99999999

		"""
		# performs initial round of two-way tracking
		closest_no_conflict = self.two_way_track()

		# performs translation correction
		self.beads_final = self.translation_correction(closest_no_conflict)

		# performs second round of two-way tracking
		closest_no_conflict = self.two_way_track()

		return closest_no_conflict

	def two_way_track(self):
		"""Performs a round of two-way tracking. Maps beads in initial state to beads in final state, then maps beads
		in final state to beads in initial state, then rectifies the two mappings to only include correspondences identified
		during both rounds of matching.

		Returns
		----------
		closest_no_conflict : np.array
			One dimensional NumPy array. Each element closest_no_conflict[i] contains the index of the bead in beads_final
			correlating to bead i in beads_init. If no corresponding bead was found, closest_no_conflict[i] = 99999999

		"""

		# forms a progress bar
		if self.print_progress:
			pbar = tqdm(total=11, desc='Tracking round', bar_format='{l_bar}{bar}|[{elapsed}<{remaining}, {rate_fmt}{postfix}]', ascii=True)
		else:
			pbar = None

		# runs the tracking
		closest_no_conflict, _ = tracking.two_way_track(self.num_feat, self.num_nearest, self.beads_init, self.beads_final, pbar=pbar) # CORRECT in tracking.py
		return closest_no_conflict

	def translation_correction(self, closest_no_conflict):
		"""Corrects for microscope translation using the method outlined in the technical overview paper. Utilizes a multivariate
		adaptive regression spline method from PyEarth to bring the average bead displacement of each z-slice to zero

		Parameters
		----------
		closest_no_conflict : np.array
			One dimensional NumPy array. Each element closest_no_conflict[i] contains the index of the bead in beads_final
			correlating to bead i in beads_init. If no corresponding bead was found, closest_no_conflict[i] = 99999999

		Returns
		----------
		beads_final_corrected : fmtrack.FMBeads
			Final bead positions after microscope translation has been corrected for
			
		"""

		self.mars_model = translation_corrector.TranslationCorrector()
		self.mars_model.create_model(self.beads_init, self.beads_final, self.cell_init, self.cell_final, self.buffer_cell, closest_no_conflict)
		beads_final_corrected = self.mars_model.correct_beads(self.beads_final)
		if self.cell_final is not None:
			self.cell_final_new = self.mars_model.correct_mesh(self.cell_final)

		return beads_final_corrected


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
