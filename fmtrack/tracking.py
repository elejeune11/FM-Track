from itertools import permutations
import matplotlib.pyplot as plt
import numpy as np
import os
from pyearth import Earth
import sys
from tqdm import tqdm
import fmtrack
np.warnings.filterwarnings('ignore') # this suppresses a FutureWarning in the PyEarth package

##########################################################################################
# import bead centers 
##########################################################################################

def import_data(file_prefix_1, file_prefix_2, root_directory):
	"""Imports cell and bead data saved using FM-Track's native file structure. Occurs after pre_processing
	step has completed

	Parameters
	----------
	file_prefix_1 : str
		File prefix for the initial state
	file_prefix_2 : str
		File prefix for the final state
	root_directory : str
		Path to folder containing folder /Gel_cell_coords

	Returns
	----------
	x_pos : np.array
		1D NumPy array specifying x positions of beads in the initial configuration
	y_pos : np.array
		1D NumPy array specifying y positions of beads in the initial configuration
	z_pos : np.array
		1D NumPy array specifying z positions of beads in the initial configuration
	x_pos_new : np.array
		1D NumPy array specifying x positions of beads in the final configuration
	y_pos_new : np.array
		1D NumPy array specifying y positions of beads in the final configuration
	z_pos_new : np.array
		1D NumPy array specifying z positions of beads in the final configuration
	cell_mesh : np.array
		NumPy array of shape (num_vertices, 3) specifying the (x,y,z) coordinates of every
		vertex in the initial cell boundary
	cell_mesh_2 : np.array
		NumPy array of shape (num_vertices, 3) specifying the (x,y,z) coordinates of every
		vertex in the final cell boundary

	"""

	cell_mesh_fname1 = root_directory + '/Gel_cell_coords/' + file_prefix_1 + '_cell_mesh.txt'
	cell_mesh_fname2 = root_directory + '/Gel_cell_coords/' + file_prefix_2 + '_cell_mesh.txt'
	cell_mesh = np.loadtxt(cell_mesh_fname1)
	cell_mesh_2 = np.loadtxt(cell_mesh_fname2)
	beads_fname1 = root_directory + '/Gel_bead_center_coords/' + file_prefix_1 + '_beads.txt'
	beads_fname2 = root_directory + '/Gel_bead_center_coords/' + file_prefix_2 + '_beads.txt'
	data1 = np.loadtxt(beads_fname1)
	data2 = np.loadtxt(beads_fname2)
	x_pos = []; y_pos = []; z_pos = [] 
	for kk in range(0,data1.shape[0]):
		x_pos.append(data1[kk,0])
		y_pos.append(data1[kk,1])
		z_pos.append(data1[kk,2])
	
	x_pos_new = []; y_pos_new = []; z_pos_new = [] 
	for kk in range(0,data2.shape[0]):
		x_pos_new.append(data2[kk,0])
		y_pos_new.append(data2[kk,1])
		z_pos_new.append(data2[kk,2])
	
	return x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, cell_mesh, cell_mesh_2

def import_data_no_cell(file_prefix_1,file_prefix_2, root_directory):
	"""Imports bead data saved using FM-Track's native file structure. Occurs after pre_processing
	step has completed. Only for trials without cell data

	Parameters
	----------
	file_prefix_1 : str
		File prefix for the initial state
	file_prefix_2 : str
		File prefix for the final state
	root_directory : str
		Path to folder containing folder /Gel_cell_coords

	Returns
	----------
	x_pos : np.array
		1D NumPy array specifying x positions of beads in the initial configuration
	y_pos : np.array
		1D NumPy array specifying y positions of beads in the initial configuration
	z_pos : np.array
		1D NumPy array specifying z positions of beads in the initial configuration
	x_pos_new : np.array
		1D NumPy array specifying x positions of beads in the final configuration
	y_pos_new : np.array
		1D NumPy array specifying y positions of beads in the final configuration
	z_pos_new : np.array
		1D NumPy array specifying z positions of beads in the final configuration

	"""

	beads_fname1 = root_directory + '/Gel_bead_center_coords/' + file_prefix_1 + '_beads.txt'
	beads_fname2 = root_directory + '/Gel_bead_center_coords/' + file_prefix_2 + '_beads.txt'
	data1 = np.loadtxt(beads_fname1)
	data2 = np.loadtxt(beads_fname2)
	x_pos = []; y_pos = []; z_pos = [] 
	for kk in range(0,data1.shape[0]):
		x_pos.append(data1[kk,0])
		y_pos.append(data1[kk,1])
		z_pos.append(data1[kk,2])
	
	x_pos_new = []; y_pos_new = []; z_pos_new = [] 
	for kk in range(0,data2.shape[0]):
		x_pos_new.append(data2[kk,0])
		y_pos_new.append(data2[kk,1])
		z_pos_new.append(data2[kk,2])
	
	return x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new 
	
##########################################################################################
# save results   
##########################################################################################

def save_res(destination, beads_init, beads_final, label_uncorrected):
	"""Saves bead positions in initial state along with displacements as text files

	Parameters
	----------
	destination : str
		Path to folder where data will be saved
	beads_init : FMBeads
		Beads in the initial position
	beads_final : FMBeads
		Beads in the final position
	label_uncorrected : bool
		Determines whether bead displacements are labeled as uncorrected or not. Intended to refer to
		whether tranlsation correction has occured

	"""

	x_pos, y_pos, z_pos = beads_init.get_xyz()
	x_pos_new, y_pos_new, z_pos_new = beads_final.get_xyz()

	if not os.path.exists(destination):
		os.makedirs(destination)

	U = x_pos_new - x_pos
	V = y_pos_new - y_pos
	W = z_pos_new - z_pos


	np.savetxt(os.path.join(destination,'X.txt'),x_pos)
	np.savetxt(os.path.join(destination,'Y.txt'),y_pos)
	np.savetxt(os.path.join(destination,'Z.txt'),z_pos)

	if label_uncorrected:
		np.savetxt(os.path.join(destination,'U_uncorrected.txt'),U)
		np.savetxt(os.path.join(destination,'V_uncorrected.txt'),V)
		np.savetxt(os.path.join(destination,'W_uncorrected.txt'),W)
	else:
		np.savetxt(os.path.join(destination,'U.txt'),U)
		np.savetxt(os.path.join(destination,'V.txt'),V)
		np.savetxt(os.path.join(destination,'W.txt'),W)

def load_res(destination,label_uncorrected):
	"""Loads bead data saved using FM-Track's native file structure (specified in save_res)

	Parameters
	----------
	destination : str
		Path to folder where data will be saved
	label_uncorrected : bool
		Determines whether bead displacements are labeled as uncorrected or not. Intended to refer to
		whether tranlsation correction has occured

	Returns
	----------
	beads_init : FMBeads
		Beads in the initial position
	beads_final : FMBeads
		Beads in the final position

	"""
	
	X = np.loadtxt(os.path.join(destination,'X.txt'))
	Y = np.loadtxt(os.path.join(destination,'Y.txt'))
	Z = np.loadtxt(os.path.join(destination,'Z.txt'))

	if label_uncorrected:
		U = np.loadtxt(os.path.join(destination,'U_uncorrected.txt'))
		V = np.loadtxt(os.path.join(destination,'V_uncorrected.txt'))
		W = np.loadtxt(os.path.join(destination,'W_uncorrected.txt'))
	else:
		U = np.loadtxt(os.path.join(destination,'U.txt'))
		V = np.loadtxt(os.path.join(destination,'V.txt'))
		W = np.loadtxt(os.path.join(destination,'W.txt'))

	beads_init_new = fmbeads.FMBeads()
	beads_final_new = fmbeads.FMBeads()

	beads_init_new.load_from_positions(X,Y,Z)
	beads_final_new.load_from_positions(X+U,Y+V,Z+W)

	return beads_init_new, beads_final_new

def get_corrected_beads(beads_init, beads_final, closest_no_conflict):
	"""Creates a FMBeads object from image data

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	x_pos : np.array
		1D NumPy array specifying x positions of beads in the initial configuration
	y_pos : np.array
		1D NumPy array specifying y positions of beads in the initial configuration
	z_pos : np.array
		1D NumPy array specifying z positions of beads in the initial configuration
	x_pos_new : np.array
		1D NumPy array specifying x positions of beads in the final configuration
	y_pos_new : np.array
		1D NumPy array specifying y positions of beads in the final configuration
	z_pos_new : np.array
		1D NumPy array specifying z positions of beads in the final configuration
	closest_no_conflict : np.array


	Returns
	----------
	X
	Y
	Z
	U
	V
	W

	"""

	x_pos, y_pos, z_pos = beads_init.get_xyz()
	x_pos_new, y_pos_new, z_pos_new = beads_final.get_xyz()

	X_init = []; Y_init = []; Z_init = []; X_final = []; Y_final = []; Z_final = [] 
	num_pts = len(x_pos)
	for kk in range(0,num_pts):
		idx = closest_no_conflict[kk]
		if idx <= num_pts:
			X_init.append(x_pos[kk])
			Y_init.append(y_pos[kk])
			Z_init.append(z_pos[kk])

			X_final.append(x_pos_new[idx])
			Y_final.append(y_pos_new[idx])
			Z_final.append(z_pos_new[idx])

	X_init = np.asarray(X_init)
	Y_init = np.asarray(Y_init)
	Z_init = np.asarray(Z_init)
	X_final = np.asarray(X_final)
	Y_final = np.asarray(Y_final)
	Z_final = np.asarray(Z_final)

	beads_init_new = fmtrack.FMBeads()
	beads_init_new.load_from_positions(X_init, Y_init, Z_init)
	beads_final_new = fmtrack.FMBeads()
	beads_final_new.load_from_positions(X_final, Y_final, Z_final)

	return beads_init_new, beads_final_new

def save_corrected_cell(destination,cell_mesh):
	np.savetxt(os.path.join(destination,'cell_mesh_2_corrected.txt'),cell_mesh.points)
	return 

##########################################################################################
# helper functions  
##########################################################################################

def get_dist_between_init_final_mat(x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, pbar=None):
	"""Calculates the distance matrix, which contains the values quantifying distance between every
	point in the initial configuration and every point in the final configuration

	Parameters
	----------
	x_pos : np.array
		1D NumPy array specifying x positions of beads in the initial configuration
	y_pos : np.array
		1D NumPy array specifying y positions of beads in the initial configuration
	z_pos : np.array
		1D NumPy array specifying z positions of beads in the initial configuration
	x_pos_new : np.array
		1D NumPy array specifying x positions of beads in the final configuration
	y_pos_new : np.array
		1D NumPy array specifying y positions of beads in the final configuration
	z_pos_new : np.array
		1D NumPy array specifying z positions of beads in the final configuration
	pbar : tqdm.tqdm, optional
		Progress bar object. Primarily for internal use when all functions are strung together.
		Recommended to not pass anything into this variable

	Returns
	----------
	dist_mat : np.array
		NumPy array of shape (num_orig, num_final) where each number is the number of beads
		in the initial and final configurations. dist_mat[i,j] specifies magnitude of vector
		pointing from point i in initial configuration to point j in final configuration.

	"""

	num_orig = len(x_pos); num_final = len(x_pos_new)
	dist_mat = np.zeros((num_orig, num_final))
	for kk in range(0,num_orig):
		x_pt = x_pos[kk]; y_pt = y_pos[kk]; z_pt = z_pos[kk] 
		dist_vec = ((x_pt - x_pos_new)**2.0 + (y_pt - y_pos_new)**2.0 + (z_pt - z_pos_new)**2.0)**(0.5)
		dist_mat[kk,:] = dist_vec
		if pbar is not None:
			pbar.update(1/num_orig)
	return dist_mat

def get_nearest_pts(num_nearest,dist_mat,pbar=None):
	"""Returns the indices of points in the final configuration considered to be "neighbors" to points in 
	the initial configuration

	Parameters
	----------
	num_nearest : int
		The number of beads to consider neighbors
	dist_mat : np.array
		NumPy array of shape (num_orig, num_final) where each number is the number of beads
		in the initial and final configurations. dist_mat[i,j] specifies magnitude of vector
		pointing from point i in initial configuration to point j in final configuration.
	pbar : tqdm.tqdm, optional
		Progress bar object. Primarily for internal use when all functions are strung together.
		Recommended to not pass anything into this variable

	Returns
	----------
	nearest_mat : np.array
		NumPy array of shape (num_points_orig, num_nearest). Each row corresponds to a specific
		bead in the initial configuration. Each row contains the indices of the closest num_nearest 
		beads in the final configuration

	"""

	num_pts_orig = dist_mat.shape[0]
	nearest_mat = np.zeros((num_pts_orig,num_nearest))
	for kk in range(0,num_pts_orig):
		args = np.argsort(dist_mat[kk,:])
		for jj in range(0,num_nearest):
			nearest_mat[kk,jj] = args[jj]
		if pbar is not None:
			pbar.update(1/num_pts_orig)
	return nearest_mat 

def get_dist_same_config(x_all,y_all,z_all,pbar=None):
	"""Get distance between points in the same configuration

	Parameters
	----------
	x_all : np.array
		1D NumPy array specifying x positions of an arbitrary set of beads
	y_all : np.array
		1D NumPy array specifying x positions of an arbitrary set of beads
	z_all : np.array
		1D NumPy array specifying x positions of an arbitrary set of beads
	pbar : tqdm.tqdm, optional
		Progress bar object. Primarily for internal use when all functions are strung together.
		Recommended to not pass anything into this variable

	Returns
	----------
	dist_mat : np.array
		NumPy array of shape (x_all.shape[0], x_all.shape[0]). dist_mat[i,j] specifies magnitude
		of vector pointing from point i to point j.

	"""

	num_pts = len(x_all)
	dist_mat = np.zeros((num_pts,num_pts))
	for kk in range(0,num_pts):
		for jj in range(kk,num_pts):
			if kk == jj:
				dist_mat[kk,jj] = 9999
			else:
				dist = ((x_all[kk] - x_all[jj])**2.0 + (y_all[kk] - y_all[jj])**2.0 + (z_all[kk] - z_all[jj])**2.0)**0.5
				dist_mat[kk,jj] = dist
				dist_mat[jj,kk] = dist
		if pbar is not None:
			pbar.update(1/num_pts)
	return dist_mat 
	
def get_feature_vectors(x_all,y_all,z_all,num_feat,pbar=None):
	"""Computes feature vector array

	Parameters
	----------
	x_all : np.array
		1D NumPy array specifying x positions of an arbitrary set of beads
	y_all : np.array
		1D NumPy array specifying x positions of an arbitrary set of beads
	z_all : np.array
		1D NumPy array specifying x positions of an arbitrary set of beads
	num_feat : int
		Number of feature vectors to compute per point
	pbar : tqdm.tqdm, optional
		Progress bar object. Primarily for internal use when all functions are strung together.
		Recommended to not pass anything into this variable

	Returns
	----------
	feat_array : np.array
		NumPy array of shape (num_pts, num_feat, 3). Each plane (feat_array[i,:,:]) represents
		a specific bead. Each row of this plane (feat_array[i,j,:]) represents a vector from 
		that bead to a neighbor (a feature vector) with x, y, and z coordinates

	"""

	num_pts = len(x_all)
	# --> get the distance between vector points 
	dist_mat = get_dist_same_config(x_all,y_all,z_all,pbar=pbar)
	# --> create the feature array
	feat_array = np.zeros((num_pts,num_feat,3))
	for kk in range(0,num_pts):
		# --> get the (#num_feat) closest points
		args = np.argsort(dist_mat[kk,:])
		for jj in range(0,num_feat):
			# --> create the feature vector array 
			feat_array[kk,jj,0] = x_all[args[jj]] - x_all[kk]
			feat_array[kk,jj,1] = y_all[args[jj]] - y_all[kk]
			feat_array[kk,jj,2] = z_all[args[jj]] - z_all[kk]
		if pbar is not None:
			pbar.update(1/num_pts)
	return feat_array 

def compute_score(feat_1,feat_2): #rows are each feature, 3 columns are x,y,z
	"""Computes the score of two feature vector sets by comparing all permutations of the two

	Parameters
	----------
	feat_1 : np.array
		NumPy array of shape (num_feat, 3) specifying a set of feature vectors for a
		particular bead
	feat_2 : np.array
		NumPy array of shape (num_feat, 3) specifying a set of feature vectors for a
		particular bead

	Returns
	----------
	best_score : float
		The score of the two feature vector sets

	"""

	num_1 = feat_1.shape[0]; num_2 = feat_2.shape[0]
	idx_1 = [kk for kk in range(0,num_1)]
	idx_2 = [kk for kk in range(0,num_2)]
	perm_list = [] 
	for p in permutations(idx_2):
		perm_list.append(p)
	score_list = [] 
	for kk in range(0,len(perm_list)):
		dist = 0
		for jj in range(0,len(perm_list[kk])):
			vec1 = feat_1[idx_1[jj],:]
			vec2 = feat_2[idx_2[perm_list[kk][jj]],:]
			dist += ((vec1[0]-vec2[0])**2.0 + (vec1[1]-vec2[1])**2.0 + (vec1[2]-vec2[2])**2.0)**0.5
		score_list.append(dist)
	best_score = np.min(score_list)
	return best_score

def get_score_info(feat_kk_orig,feat_array_fin,feat_idx_fin):
	"""Matching scores

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	feat_kk_orig : np.array
		Original beads feature array of shape (num_pts, num_feat, 3). Each plane (feat_array[i,:,:]) 
		represents a specific bead. Each row of this plane (feat_array[i,j,:]) represents a vector from 
		that bead to a neighbor (a feature vector) with x, y, and z coordinates
	feat_array_fin : np.array
		Final beads feature array
	feat_idx_fin : np.array
		1D NumPy array specifying the indices that order beads_final according to how they correlate with 
		beads_initial

	Returns
	----------
	idx_all : 
	score_all

	"""

	num_nearest = feat_idx_fin.shape[0] 
	vec = np.zeros((num_nearest))
	for kk in range(0,num_nearest):
		# --> get score
		feat_1 = feat_kk_orig
		feat_2 = feat_array_fin[int(feat_idx_fin[kk]),:,:]
		dist_score = compute_score(feat_1,feat_2)
		vec[kk] = dist_score
	idx_sort = np.argsort(vec)
	idx_all = np.zeros((num_nearest))
	score_all = np.zeros((num_nearest))
	for jj in range(0,num_nearest):
		idx_all[jj] = feat_idx_fin[idx_sort[jj]]
		score_all[jj] = vec[idx_sort[jj]]
	return idx_all, score_all 

def get_closest_features(feat_array_orig, feat_array_fin, nearest_mat, num_feat,num_nearest,pbar=None):	
	"""Comparing feature vectors

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	feat_array_orig
	feat_array_fin
	nearest_mat
	num_feat
	num_nearest
	pbar

	Returns
	----------
	matching_idx_sorted
	matching_score_sorted

	"""

	# --> get the matching scores 
	num_pts = feat_array_orig.shape[0]  
	matching_idx_sorted = np.zeros((num_pts,num_nearest))
	matching_score_sorted = np.zeros((num_pts,num_nearest))

	for kk in range(0,num_pts):
		feat_kk_orig = feat_array_orig[kk,:,:]
		feat_idx_fin = nearest_mat[kk,:] 
		idx_all, score_all = get_score_info(feat_kk_orig,feat_array_fin,feat_idx_fin)
		matching_idx_sorted[kk,:] = idx_all
		matching_score_sorted[kk,:] = score_all 
		if pbar is not None:
			pbar.update(1/num_pts)
	
	return matching_idx_sorted, matching_score_sorted

def iterate_closest(matching_idx_sorted, matching_score_sorted, num_nearest, pbar=None):
	"""Iterate through  

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	matching_idx_sorted
	matching_score_sorted
	num_nearest

	Returns
	----------
	closest_no_conflict

	"""

	closest = matching_idx_sorted[:,0]
	num_beads = matching_idx_sorted.shape[0] 
	idx_resolved = np.zeros((matching_idx_sorted.shape[0]))
	for zz in range(0,num_nearest - 1):
		conflict_counter = 0 
		closest_no_conflict = np.zeros((closest.shape[0]))
		for kk in range(0,num_beads):
			# --> check to see if there is a conflict
			condition = closest == kk
			vals = np.nonzero(condition)[0]
			if vals.shape[0] > 1: # <-- indicates that there is a conflict 
				conflict_counter += 1 
				# identify which of the conflicts has the highest score, keep this one, move the rest to the next-best match
				scores = [] 
				for zz in range(0,vals.shape[0]):
					scores.append(matching_score_sorted[vals[zz],int(idx_resolved[vals[zz]])])
				args = np.argsort(scores)
				# keep the one that has the highest score
				closest_no_conflict[vals[args[0]]] = kk
				# re-assign the rest to the next highest score, update indexing
				for zz in range(1,vals.shape[0]):
					closest_no_conflict[vals[args[zz]]] = matching_score_sorted[ vals[args[zz]], int(idx_resolved[vals[args[zz]]])+1 ]
					idx_resolved[vals[args[zz]]] = idx_resolved[vals[args[zz]]] + 1 
			elif vals.shape[0] == 1:
				closest_no_conflict[int(vals[0])] = kk 
			if pbar is not None:
				pbar.update(1/(num_nearest-1)/num_beads)
		closest = np.copy(closest_no_conflict) 
	return closest_no_conflict 

def check_backward_forward(closest_no_conflict_f,closest_no_conflict_b):
	"""Check the results of running the tracking algorithm both backwards and forwards, discard results that don't agree

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	closest_no_conflict_f
	closest_no_conflict_b

	Returns
	----------
	closest_no_conflict
	idx_ignored

	"""

	num_pts = len(closest_no_conflict_f)
	idx_ignored = [] 
	closest_no_conflict = [] 

	for kk in range(0,num_pts):
		match_f = int(closest_no_conflict_f[kk])
		match_b = int(closest_no_conflict_b[match_f])
		if match_b == kk:
			# good match 
			closest_no_conflict.append(match_f)
		else:
			# bad match 
			closest_no_conflict.append(99999999)
			idx_ignored.append(kk)
	return closest_no_conflict, idx_ignored


def translation_correction(cell_init, cell_final, buffer_cell,\
	x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict):
	"""Sometimes, microscope translation and/or gel swelling can lead to translation
	This function will correct for this translation and then re-run the tracking algorithm
	This only deals with corrections with respect to the z direction 

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	cell_init
	cell_final
	buffer_cell
	x_pos
	y_pos
	z_pos
	x_pos_new
	y_pos_new
	z_pos_new
	closest_no_conflict

	Returns
	----------
	x_pos_new
	y_pos_new
	z_pos_new
	cell_final_new
	fig

	"""

	points_init = cell_init.points
	points_final = cell_final.points
	
	x_min = np.min([np.min(points_init[:,0]),np.min(points_final[:,0])]) - buffer_cell 
	x_max = np.max([np.max(points_init[:,0]),np.max(points_final[:,0])]) + buffer_cell
	y_min = np.min([np.min(points_init[:,1]),np.min(points_final[:,1])]) - buffer_cell
	y_max = np.max([np.max(points_init[:,1]),np.max(points_final[:,1])]) + buffer_cell
	z_min = np.min([np.min(points_init[:,2]),np.min(points_final[:,2])]) - buffer_cell
	z_max = np.max([np.max(points_init[:,2]),np.max(points_final[:,2])]) + buffer_cell
	
	num_pts = len(x_pos)
	X = []; Y = []; Z = []; U = []; V = []; W = [] 
	for kk in range(0,num_pts):
		idx = closest_no_conflict[kk]
		if idx < len(closest_no_conflict):
			U.append(x_pos_new[idx] - x_pos[kk])
			V.append(y_pos_new[idx] - y_pos[kk])
			W.append(z_pos_new[idx] - z_pos[kk])
			X.append(x_pos_new[idx]); Y.append(y_pos_new[idx]); Z.append(z_pos_new[idx])
	
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
	
	# --> fit MARS models 
	model_U = Earth(max_degree=2,max_terms=10)
	model_U.fit(Z_safe,U_safe)
	model_V = Earth(max_degree=2,max_terms=10)
	model_V.fit(Z_safe,V_safe)
	model_W = Earth(max_degree=2,max_terms=10)
	model_W.fit(Z_safe,W_safe)
		
	# --> re-define Z 
	pred_U = model_U.predict(z_pos_new)
	pred_V = model_V.predict(z_pos_new)
	pred_W = model_W.predict(z_pos_new)
	
	# --> correct new bead positions 
	for kk in range(0,len(x_pos_new)):
		x_pos_new[kk] = x_pos_new[kk] - pred_U[kk] 
		y_pos_new[kk] = y_pos_new[kk] - pred_V[kk]
		z_pos_new[kk] = z_pos_new[kk] - pred_W[kk] 
	
	# --> correct new cell position 
	pred_cell_0 = model_U.predict(points_final[:,0])
	pred_cell_1 = model_V.predict(points_final[:,1])
	pred_cell_2 = model_W.predict(points_final[:,2])
	
	points_final_new = np.zeros(points_final.shape)
	points_final_new[:,0] = points_final[:,0] - pred_cell_0
	points_final_new[:,1] = points_final[:,1] - pred_cell_1
	points_final_new[:,2] = points_final[:,2] - pred_cell_2

	cell_final_new = cell_final
	cell_final_new.points = points_final_new
	
	# --> plot MARS models 
	Z_line = np.linspace(np.min(Z),np.max(Z),100)
	pred_line_U = model_U.predict(Z_line)
	pred_line_V = model_V.predict(Z_line)
	pred_line_W = model_W.predict(Z_line)

	fig = plt.figure(figsize=(15,5))
	fig.tight_layout()
	ax1 = fig.add_subplot(1,3,1)
	ax1.plot(Z,U,'b.',label='x raw')
	ax1.plot(Z_line,pred_line_U,'k--',label='fit')
	ax1.set_xlabel('z position'); ax1.set_ylabel('displacement')
	ax1.legend(); ax1.set_title('x displacements')
	ax2 = fig.add_subplot(1,3,2)
	ax2.plot(Z,V,'r.',label='y raw')
	ax2.plot(Z_line,pred_line_V,'k--',label='fit')
	ax2.set_xlabel('z position'); ax2.set_ylabel('displacement')
	ax2.legend(); ax2.set_title('y displacements')
	ax3 = fig.add_subplot(1,3,3)
	ax3.plot(Z,W,'g.',label='z raw')
	ax3.plot(Z_line,pred_line_W,'k--',label='fit')
	ax3.set_xlabel('z position'); ax3.set_ylabel('displacement')
	ax3.legend(); ax3.set_title('z displacements')

	return x_pos_new, y_pos_new, z_pos_new, cell_final_new, fig

##########################################################################################
# main tracking functions 
##########################################################################################

def two_way_track(num_feat, num_nearest, beads_init, beads_final, pbar=None): # (coded to not re-compute distances when possible )
	"""Two way tracking 

	#################
	NEEDS DESCRIPTION
	#################

	Parameters
	----------
	num_feat
	num_nearest
	x_pos
	y_pos
	z_pos
	x_pos_new
	y_pos_new
	z_pos_new
	pbar

	Returns
	----------
	closest_no_conflict
	idx_ignored

	"""

	x_pos, y_pos, z_pos = beads_init.get_xyz()
	x_pos_new, y_pos_new, z_pos_new = beads_final.get_xyz()

	# --> get the distance between each point initial to final config set up 
	dist_mat_i_to_f = get_dist_between_init_final_mat(x_pos,y_pos,z_pos,x_pos_new,y_pos_new,z_pos_new, pbar=pbar) 
	dist_mat_f_to_i = dist_mat_i_to_f.T 
	
	# --> get the feature vectors for both configurations
	feat_array_i = get_feature_vectors(x_pos,y_pos,z_pos,num_feat, pbar=pbar)
	feat_array_f = get_feature_vectors(x_pos_new,y_pos_new,z_pos_new,num_feat, pbar=pbar)
	
	# --> get the nearest matrices 
	nearest_mat_i_to_f = get_nearest_pts(num_nearest,dist_mat_i_to_f, pbar=pbar)
	nearest_mat_f_to_i = get_nearest_pts(num_nearest,dist_mat_f_to_i, pbar=pbar)
	
	# --> matching indicies and score forward
	matching_idx_sorted_f, matching_score_sorted_f = get_closest_features(feat_array_i, feat_array_f, nearest_mat_i_to_f, num_feat, num_nearest, pbar=pbar)
	
	# --> matching indicies and score backward 
	matching_idx_sorted_b, matching_score_sorted_b = get_closest_features(feat_array_f, feat_array_i, nearest_mat_f_to_i, num_feat, num_nearest, pbar=pbar)

	# --> forward run
	closest_no_conflict_f = iterate_closest(matching_idx_sorted_f, matching_score_sorted_f, num_nearest, pbar=pbar)
	
	# --> backward run
	closest_no_conflict_b = iterate_closest(matching_idx_sorted_b, matching_score_sorted_b, num_nearest, pbar=pbar)

	if pbar is not None:
		#pbar.refresh()
		pbar.close()
	
	# --> return matched
	closest_no_conflict, idx_ignored = check_backward_forward(closest_no_conflict_f,closest_no_conflict_b)
	
	return closest_no_conflict, idx_ignored