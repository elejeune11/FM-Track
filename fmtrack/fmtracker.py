import numpy as np 
import os
from . import tracking
from . import fmplot
from . import post_process
from . import fmbeads
from . import fmmesh
import fmtrack
from . import translation_corrector
import pickle
from tqdm import tqdm
from scipy import stats

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
		self.beads_final_uncorrected = None
		self.mars_model = None

		self.closest_no_conflict = None

		self.print_progress = True

		self.use_box = True
		self.num_feat = 5
		self.num_nearest = 15
		self.buffer_cell = 30
		self.track_type = 2 # type 1 will NOT perform translation correction, type 2 will

	def save(self,filename):
		# saves entire FMTracker object by pickling it
		pickle.dump(self, open(filename,'wb'))

	def run_tracking(self):
		"""Main tracking function. If self.track_type is 2, run_tracking() runs the more robust tracking algorithm,
		correcting for microscope translation. If self.track_type is 1, run_tracking() does not correct for microscope
		translation.

		"""
			
		if self.track_type == 1:
			closest_no_conflict = self.two_way_track()
		elif self.track_type == 2:
			closest_no_conflict = self.track_correct_track()

		self.beads_init_new, self.beads_final_new, self.matches = tracking.get_corrected_beads(self.beads_init, self.beads_final, closest_no_conflict, return_matches=True)

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
		self.beads_final_uncorrected = fmtrack.FMBeads()
		self.beads_final_uncorrected.points = np.copy(self.beads_final.points)
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
		closest_no_conflict, _ = tracking.two_way_track(self.num_feat, self.num_nearest, self.beads_init, self.beads_final, pbar=pbar)
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
		self.mars_model.create_model(self.beads_init, self.beads_final, self.cell_init, self.cell_final, self.buffer_cell, closest_no_conflict, use_box=self.use_box)
		beads_final_corrected = self.mars_model.correct_beads(self.beads_final)
		if self.cell_final is not None:
			self.cell_final_new = self.mars_model.correct_mesh(self.cell_final)

		return beads_final_corrected

	def remove_spurious(self, z_score_spurious=2, return_isspurious=False):
		"""Removes spurious bead matches. "Spurious" beads fulfill these two requirements:
		
		(1) They are in the "far field," meaning they are further from the cell boundary
			than some distance threshold self.buffer_cell
		(2) Of all the far-field bead displacements, they deviate from mean displacement
			by z_score_spurious standard deviations

        """

		# identify far field beads
		is_farfield = identify_farfield(self.beads_init_new.points, self.cell_init.points, self.buffer_cell)
		idx_farfield = np.argwhere(is_farfield).flatten()
		disps_all = self.beads_final_new.points - self.beads_init_new.points
		disps_farfield = self.beads_final_new.points[idx_farfield,:] - self.beads_init_new.points[idx_farfield,:] 

		# compute spurious displacements in u, v, and w directions
		is_spurious = np.zeros(is_farfield.shape, dtype=bool)
		for i in range(3):
			mean = np.mean(disps_farfield[:,i])
			std = np.std(disps_farfield[:,i])
			z_score_array = (disps_all[:,i] - mean) / std
			spurious = np.logical_or(z_score_array < -z_score_spurious, z_score_array > z_score_spurious)
			is_spurious = np.logical_or(is_spurious, spurious)

		# must be both spurious and farfield for consideration
		is_spurious = np.logical_and(is_spurious, is_farfield)
		idx_spurious = np.argwhere(is_spurious).flatten()
		idx_nonspurious = np.argwhere(np.logical_not(is_spurious)).flatten()

		# return the non-spurious beads
		self.beads_init_new.points = self.beads_init_new.points[idx_nonspurious,:]
		self.beads_final_new.points = self.beads_final_new.points[idx_nonspurious,:]

		if return_isspurious:
			return is_spurious

	def save_all(self, folderpath):
		# saves the meshes, bead positions, and relevant graphs
		os.makedirs(folderpath, exist_ok=True)

		# saves the tracker itself
		self.save(folderpath+'/tracker.sav')

		# saves the cells as msh files
		if self.cell_init is not None:
			self.cell_init.save_msh_file(folderpath+'/cell_init.msh')
		if self.cell_final is not None:
			self.cell_final.save_msh_file(folderpath+'/cell_final.msh')

		# saves all of the bead positions, even those not in a match
		path = folderpath+'/all_beads'
		os.makedirs(path, exist_ok=True)
		np.savetxt(path+'/beads_init.txt', self.beads_init.points)
		if self.beads_final_uncorrected is None:
			np.savetxt(path+'/beads_final.txt', self.beads_final.points)
		else:
			np.savetxt(path+'/beads_final_uncorrected.txt', self.beads_final_uncorrected.points)
			np.savetxt(path+'/beads_final_corrected.txt', self.beads_final.points)

		# saves matched bead positions
		path = folderpath+'/matched_beads'
		os.makedirs(path, exist_ok=True)
		np.savetxt(path+'/beads_init.txt', self.beads_init_new.points)
		np.savetxt(path+'/beads_final.txt', self.beads_final_new.points)
		np.savetxt(path+'/matches.txt', self.matches)

		# saves graphs
		plotter = fmtrack.FMPlot(self)
		plotter.save_native_plots(folderpath)


	def save_res(self,folder,label_uncorrected):
		# saves bead positions in initial state along with displacements as text files
		tracking.save_res(folder, self.beads_init_new, self.beads_final_new, label_uncorrected)

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
'''
def identify_nearfield(points, cell, thresh):
	idx_nearfield = np.array([])
	for i in range(points.shape[0]):
		if np.any(np.linalg.norm(cell - points[i,:], axis=1) < thresh):
			idx_nearfield = np.append(idx_farfield, i)

	return idx_nearfield.astype(int)
'''
def identify_farfield(points, cell, thresh):
	idx_nearfield = np.zeros(points.shape[0]).astype(bool)
	for i in range(points.shape[0]):
		idx_nearfield[i] = np.any(np.linalg.norm(cell - points[i,:], axis=1) < thresh)

	return np.logical_not(idx_nearfield)