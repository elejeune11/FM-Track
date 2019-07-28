from itertools import permutations
import matplotlib.pyplot as plt
import numpy as np
import os
from pyearth import Earth
import sys
np.warnings.filterwarnings('ignore') # this suppresses a FutureWarning in the PyEarth package
##########################################################################################
# import bead centers 
##########################################################################################
def import_data(file_prefix_1, file_prefix_2):
	cell_mesh_fname1 = 'Gel_cell_coords/' + file_prefix_1 + '_cell_mesh.txt'
	cell_mesh_fname2 = 'Gel_cell_coords/' + file_prefix_2 + '_cell_mesh.txt'
	cell_mesh = np.loadtxt(cell_mesh_fname1)
	cell_mesh_2 = np.loadtxt(cell_mesh_fname2)
	beads_fname1 = 'Gel_bead_center_coords/' + file_prefix_1 + '_beads.txt'
	beads_fname2 = 'Gel_bead_center_coords/' + file_prefix_2 + '_beads.txt'
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

def import_data_no_cell(file_prefix_1,file_prefix_2):
	beads_fname1 = 'Gel_bead_center_coords/' + file_prefix_1 + '_beads.txt'
	beads_fname2 = 'Gel_bead_center_coords/' + file_prefix_2 + '_beads.txt'
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
def save_res(destination, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict, label_uncorrected):
	X = []; Y = []; Z = []; U = []; V = []; W = [] 
	num_pts = len(x_pos)
	for kk in range(0,num_pts):
		idx = closest_no_conflict[kk]
		if idx <= num_pts:
			X.append(x_pos[kk])
			Y.append(y_pos[kk])
			Z.append(z_pos[kk])
			U.append(x_pos_new[idx] - x_pos[kk])
			V.append(y_pos_new[idx] - y_pos[kk])
			W.append(z_pos_new[idx] - z_pos[kk]) 

	np.savetxt(destination + '/X.txt',np.asarray(X))
	np.savetxt(destination + '/Y.txt',np.asarray(Y))
	np.savetxt(destination + '/Z.txt',np.asarray(Z))

	if label_uncorrected:
		np.savetxt(destination + '/U_uncorrected.txt',np.asarray(U))
		np.savetxt(destination + '/V_uncorrected.txt',np.asarray(V))
		np.savetxt(destination + '/W_uncorrected.txt',np.asarray(W))
	else:
		np.savetxt(destination + '/U.txt',np.asarray(U))
		np.savetxt(destination + '/V.txt',np.asarray(V))
		np.savetxt(destination + '/W.txt',np.asarray(W))
	return

def save_corrected_cell(destination,cell_mesh):
	np.savetxt(destination + '/cell_mesh_2_corrected.txt',cell_mesh)
	return 
##########################################################################################
# helper functions  
##########################################################################################
# --> finding a point's neighborhood 
# every point in the initial config is compared to every point in the final config
# distance matrix is #original (rows) x #final (columns)
def get_dist_between_init_final_mat(x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new):
	num_orig = len(x_pos); num_final = len(x_pos_new)
	dist_mat = np.zeros((num_orig, num_final))
	for kk in range(0,num_orig):
		x_pt = x_pos[kk]; y_pt = y_pos[kk]; z_pt = z_pos[kk] 
		dist_vec = ((x_pt - x_pos_new)**2.0 + (y_pt - y_pos_new)**2.0 + (z_pt - z_pos_new)**2.0)**(0.5)
		dist_mat[kk,:] = dist_vec
	return dist_mat

# --> finding a point's neighborhood 
#return the nearest pt in the final configuration to a given pt in the initial config
def get_nearest_pts(num_nearest,dist_mat):
	num_pts_orig = dist_mat.shape[0]
	nearest_mat = np.zeros((num_pts_orig,num_nearest))
	for kk in range(0,num_pts_orig):
		args = np.argsort(dist_mat[kk,:])
		for jj in range(0,num_nearest):
			nearest_mat[kk,jj] = args[jj]
	return nearest_mat 

# --> get distance between points in the same configuration
def get_dist_same_config(x_all,y_all,z_all):
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
	return dist_mat 
	
# --> computing feature vectors (3D array)
def get_feature_vectors(x_all,y_all,z_all,num_feat):
	num_pts = len(x_all)
	# --> get the distance between vector points 
	dist_mat = get_dist_same_config(x_all,y_all,z_all)
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
	return feat_array 

# --> compute score 
def compute_score(feat_1,feat_2): #rows are each feature, 3 columns are x,y,z
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

# --> matching scores
def get_score_info(feat_kk_orig,feat_array_fin,feat_idx_fin):
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

# --> comparing feature vectors
def get_closest_features(feat_array_orig, feat_array_fin, nearest_mat, num_feat,num_nearest):	
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
	
	return matching_idx_sorted, matching_score_sorted


# --> iterate through  
def iterate_closest(matching_idx_sorted, matching_score_sorted, num_nearest):
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
		closest = np.copy(closest_no_conflict) 
	return closest_no_conflict 

# --> check the results of running the tracking algorithm both backwards and forwards, 
#			discard results that don't agree
def check_backward_forward(closest_no_conflict_f,closest_no_conflict_b):	
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

# --> sometimes, microscope translation and/or gel swelling can lead to translation
# --> this function will correct for this translation and then re-run the tracking algorithm
# --> this only deals with corrections with respect to the z direction 
def translation_correction(cell_mesh, cell_mesh_2, buffer_cell,\
	x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict, directory ):
	
	x_min = np.min([np.min(cell_mesh[:,0]),np.min(cell_mesh_2[:,0])]) - buffer_cell 
	x_max = np.max([np.max(cell_mesh[:,0]),np.max(cell_mesh_2[:,0])]) + buffer_cell
	y_min = np.min([np.min(cell_mesh[:,1]),np.min(cell_mesh_2[:,1])]) - buffer_cell
	y_max = np.max([np.max(cell_mesh[:,1]),np.max(cell_mesh_2[:,1])]) + buffer_cell
	z_min = np.min([np.min(cell_mesh[:,2]),np.min(cell_mesh_2[:,2])]) - buffer_cell
	z_max = np.max([np.max(cell_mesh[:,2]),np.max(cell_mesh_2[:,2])]) + buffer_cell
	
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
	pred_cell_0 = model_U.predict(cell_mesh_2[:,0])
	pred_cell_1 = model_V.predict(cell_mesh_2[:,1])
	pred_cell_2 = model_W.predict(cell_mesh_2[:,2])
	
	cell_mesh_2_new = np.zeros(cell_mesh_2.shape)
	cell_mesh_2_new[:,0] = cell_mesh_2[:,0] - pred_cell_0
	cell_mesh_2_new[:,1] = cell_mesh_2[:,1] - pred_cell_1
	cell_mesh_2_new[:,2] = cell_mesh_2[:,2] - pred_cell_2
	
	# --> plot MARS models 
	Z_line = np.linspace(np.min(Z),np.max(Z),100)
	pred_line_U = model_U.predict(Z_line)
	pred_line_V = model_V.predict(Z_line)
	pred_line_W = model_W.predict(Z_line)
	
	plt.figure(figsize=(15,5))
	plt.subplot(1,3,1)
	plt.plot(Z,U,'b.',label='x raw')
	plt.plot(Z_line,pred_line_U,'k--',label='fit')
	plt.xlabel('z position'); plt.ylabel('displacement')
	plt.tight_layout(); plt.legend(); plt.title('x displacements')
	plt.subplot(1,3,2)
	plt.plot(Z,V,'r.',label='y raw')
	plt.plot(Z_line,pred_line_V,'k--',label='fit')
	plt.xlabel('z position'); plt.ylabel('displacement')
	plt.tight_layout(); plt.legend(); plt.title('y displacements')
	plt.subplot(1,3,3)
	plt.plot(Z,W,'g.',label='z raw')
	plt.plot(Z_line,pred_line_W,'k--',label='fit')
	plt.xlabel('z position'); plt.ylabel('displacement')
	plt.tight_layout(); plt.legend(); plt.title('z displacements')
	plt.savefig(directory + '/translation_correction.png')
	
	return x_pos_new, y_pos_new, z_pos_new, cell_mesh_2_new 

##########################################################################################
# main tracking functions 
##########################################################################################
# --> two way tracking 
def two_way_track(num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new): # (coded to not re-compute distances when possible )
	# --> get the distance between each point initial to final config set up 
	dist_mat_i_to_f = get_dist_between_init_final_mat(x_pos,y_pos,z_pos,x_pos_new,y_pos_new,z_pos_new) 
	dist_mat_f_to_i = dist_mat_i_to_f.T 
	
	# --> get the feature vectors for both configurations
	feat_array_i = get_feature_vectors(x_pos,y_pos,z_pos,num_feat)
	feat_array_f = get_feature_vectors(x_pos_new,y_pos_new,z_pos_new,num_feat)
	
	# --> get the nearest matrices 
	nearest_mat_i_to_f = get_nearest_pts(num_nearest,dist_mat_i_to_f)
	nearest_mat_f_to_i = get_nearest_pts(num_nearest,dist_mat_f_to_i)
	
	# --> matching indicies and score forward
	matching_idx_sorted_f, matching_score_sorted_f = get_closest_features(feat_array_i, feat_array_f, nearest_mat_i_to_f, num_feat, num_nearest)
	
	# --> matching indicies and score backward 
	matching_idx_sorted_b, matching_score_sorted_b = get_closest_features(feat_array_f, feat_array_i, nearest_mat_f_to_i, num_feat, num_nearest)
	
	# --> forward run
	closest_no_conflict_f = iterate_closest(matching_idx_sorted_f, matching_score_sorted_f, num_nearest)
	
	# --> backward run
	closest_no_conflict_b = iterate_closest(matching_idx_sorted_b, matching_score_sorted_b, num_nearest)
	
	# --> return matched
	closest_no_conflict, idx_ignored = check_backward_forward(closest_no_conflict_f,closest_no_conflict_b)
	
	return closest_no_conflict, idx_ignored

# --> two way tracking with translation correction (track, correct, track again)
def track_correct_track(folder, num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, directory, cell_mesh, cell_mesh_2, buffer_cell):
	# --> run two way tracking 
	closest_no_conflict, idx_ignored = two_way_track(num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new)
	
	label_uncorrected = True
	save_res(folder, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict,label_uncorrected)
	# --> correct translation 	
	x_pos_new, y_pos_new, z_pos_new, cell_mesh_2_new = translation_correction(cell_mesh, cell_mesh_2, buffer_cell,\
		x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict, directory )
	
	# --> run two way tracking 
	closest_no_conflict, idx_ignored = two_way_track(num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new)
	
	return closest_no_conflict, idx_ignored, x_pos_new, y_pos_new, z_pos_new, cell_mesh_2_new 

# --> main tracking function 
def track_main_call(type,file_prefix_1, file_prefix_2, num_feat, num_nearest, buffer_cell):
	folder = 'Track_' + file_prefix_1 + '_to_' + file_prefix_2 
	# --> import files 
	x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, cell_mesh, cell_mesh_2 = \
		import_data(file_prefix_1, file_prefix_2)
	
	# --> create directory and copy files into it
	if not os.path.exists(folder):
		os.makedirs(folder)
		
	if type == 1:
		# --> two way tracking
		closest_no_conflict, idx_ignored = two_way_track(num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new)
	elif type == 2:
		# --> track, correct, track 
		closest_no_conflict, idx_ignored, x_pos_new, y_pos_new, z_pos_new, cell_mesh_2 =\
			track_correct_track(folder, num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, folder, cell_mesh, cell_mesh_2, buffer_cell)
		save_corrected_cell(folder,cell_mesh_2)
	# --> save tracking results in the proper directory 
	label_uncorrected = False
	save_res(folder, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict,label_uncorrected)
	return closest_no_conflict

# --> main tracking function, no cell, primary utility is for de-bugging and performance checks w/ synthetic data 
def track_no_cell(type,file_prefix_1,file_prefix_2,num_feat,num_nearest):
	folder = 'Track_' + file_prefix_1 + '_to_' + file_prefix_2 
	x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new = import_data_no_cell(file_prefix_1,file_prefix_2)
	
	# --> create directory and copy files into it
	if not os.path.exists(folder):
		os.makedirs(folder)
	
	if type == 1:
		# --> two way tracking
		closest_no_conflict, idx_ignored = two_way_track(num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new)
	elif type == 2:
		# --> dummy cell mesh
		cell_mesh = np.zeros((2,3))
		cell_mesh_2 = np.zeros((2,3))
		buffer_cell = 0 
		# --> track, correct, track 
		closest_no_conflict, idx_ignored, x_pos_new, y_pos_new, z_pos_new, cell_mesh_2 =\
			track_correct_track(folder, num_feat, num_nearest, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, folder, cell_mesh, cell_mesh_2, buffer_cell)
		
	# --> save tracking results in the proper directory 
	label_uncorrected = False
	save_res(folder, x_pos, y_pos, z_pos, x_pos_new, y_pos_new, z_pos_new, closest_no_conflict,label_uncorrected)
	return closest_no_conflict
	