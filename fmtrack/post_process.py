import fmtrack
import os
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (RBF, Matern, RationalQuadratic,
                                              ExpSineSquared, DotProduct,
                                              ConstantKernel, WhiteKernel)
from sklearn.neighbors import KernelDensity
from sklearn import preprocessing
##########################################################################################
# get filepath for matplotlib style
##########################################################################################
stylepath = os.path.dirname(os.path.abspath(fmtrack.__file__)) + '/el_papers.mplstyle'

##########################################################################################
# import data 
##########################################################################################
def import_cell_info(file_prefix_1,file_prefix_2,root_directory):
	cell_mesh_1 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_1 + '_cell_mesh.txt')
	cell_normal_1 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_1 + '_cell_normals.txt')
	cell_center_1 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_1 + '_cell_center.txt')
	cell_vol_1 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_1 + '_cell_volume.txt')
	cell_mesh_2 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_2 + '_cell_mesh.txt')
	cell_normal_2 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_2 + '_cell_normals.txt')
	cell_center_2 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_2 + '_cell_center.txt')
	cell_vol_2 = np.loadtxt(root_directory+'/Gel_cell_coords/' + file_prefix_2 + '_cell_volume.txt')
	return cell_mesh_1, cell_normal_1, cell_center_1, cell_vol_1, cell_mesh_2, cell_normal_2, cell_center_2, cell_vol_2
	
def import_bead_disps(folder):
	X = np.loadtxt(folder + '/X.txt')
	Y = np.loadtxt(folder + '/Y.txt')
	Z = np.loadtxt(folder + '/Z.txt')
	U = np.loadtxt(folder + '/U.txt')
	V = np.loadtxt(folder + '/V.txt')
	W = np.loadtxt(folder + '/W.txt')
	return X, Y, Z, U, V, W 
##########################################################################################
# additional computations based on the data 
##########################################################################################
# can be implemented as needed based on some rule for excluding outliers 
def remove_outliers(X, Y, Z, U, V, W):
	# maximum plausible displacement
	# z-score based displacement 
	# more complex strategy for determining outliers (hot-spot analysis)
	# this code should be implemented on a case by case basis 
	return X_new, Y_new, Z_new, U_new, V_new, W_new

# compare bead displacement to it's neighbors
def color_point_neighbor_similarity(X, Y, Z, U, V, W, num_neigh):
	num_beads = X.shape[0]
	neigh_score = [] 
	num_pts = X.shape[0] 
	for kk in range(0,num_pts):
		x = X[kk]; y = Y[kk]; z = Z[kk]
		u = U[kk]; v = V[kk]; w = W[kk]
		dist_all = ((x - X)**2.0 + (y - Y)**2.0 + (z - Z)**2.0)**(1.0/2.0)
		dist_all_sorted = np.argsort(dist_all)
		score_dist = np.zeros((num_neigh))
		for jj in range(0,num_neigh):
			idx = dist_all_sorted[jj]
			u2 = U[idx]; v2 = V[idx]; w2 = W[idx]
			score_dist[jj] = ((u - u2)**2.0 + (v - v2)**2.0 + (w - w2)**2.0)**(1.0/2.0)
		neigh_score.append(np.mean(score_dist))
	return neigh_score 

# compare bead displacement direction to the initial cell configuration 
def color_point_direction(X, Y, Z, U, V, W, cell_mesh, cell_normal):
	num_beads = X.shape[0]
	dir_score = [] 
	dist_from_cell = [] 
	mag_list = []
	# --> down sample the cell mesh (computational efficiency)
	num_pts = X.shape[0] 
	samp = np.random.randint(cell_mesh.shape[0]-1,size=np.min([num_pts,10000]))
	reduced_cell_mesh = cell_mesh[samp,:]
	reduced_cell_normal = cell_normal[samp,:]
	for kk in range(0,num_pts):
		x = X[kk]; y = Y[kk]; z = Z[kk]
		u = U[kk]; v = V[kk]; w = W[kk]
		mag = (u**2.0 + v**2.0 + w**2.0)**(1.0/2.0)
		du = u/mag; dv = v/mag; dw = w/mag
		dist_all = ((x - reduced_cell_mesh[:,0])**2.0 + (y - reduced_cell_mesh[:,1])**2.0\
			+ (z - reduced_cell_mesh[:,2])**2.0)**(1.0/2.0)
		arg = np.argmin(dist_all)
		val = du*reduced_cell_normal[arg,0] + dv*reduced_cell_normal[arg,1] + dw*reduced_cell_normal[arg,2]
		dir_score.append(val)
		dist_from_cell.append(dist_all[arg])
		mag_list.append(mag)
	return dir_score, dist_from_cell, mag_list 

# compute bead displacement to the domain edge 
def compute_dist_from_edge(X, Y, Z, X_DIM, Y_DIM, Z_DIM):
	num_pts = X.shape[0]
	dist_from_edge = [] 
	for kk in range(0,num_pts):
		x = X[kk]; y = Y[kk]; z = Z[kk]
		x_edge = np.min([np.abs(x),np.abs(X_DIM-x)])
		y_edge = np.min([np.abs(y),np.abs(Y_DIM-y)])
		z_edge = np.min([np.abs(z),np.abs(Z_DIM-z)])
		dist = np.min([x_edge,y_edge,z_edge])
		dist_from_edge.append(dist)
	return dist_from_edge

# bin data to assist with plotting 
def mean_bins(data1, data2):
	cent_mark = [10,30,50,70]
	less_than = [20,40,60,80]
	mean_val = [] 
	arg = np.argsort(data1)
	data1 = np.sort(data1)
	data2 = data2[arg]
	idx_d = 0 
	for idx_l in range(0,len(less_than)):
		arr = []; arr.append(0)
		while idx_d < data1.shape[0] and data1[idx_d] < less_than[idx_l]:
			arr.append(data2[idx_d])
			idx_d +=  1 
		mean_val.append(np.mean(arr)) 
	return cent_mark, mean_val
##########################################################################################
# plot raw data (position)
##########################################################################################
# --> y axis is bead displacement magnitude, x axis is distance from cell surface 
def plot_surface_disp(axi,cell_mesh,dist_from_edge,dist_from_cell, mag_list):
	#--> remove points 
	keep = [] 
	for kk in range(0,len(dist_from_edge)):
		if dist_from_edge[kk] > 5:
			keep.append(kk)
	keep = np.asarray(keep)
	dist_from_cell = np.asarray(dist_from_cell)
	mag_list = np.asarray(mag_list)
	cent_mark,mean_val = mean_bins(dist_from_cell[keep],mag_list[keep])
	axi.plot(dist_from_cell[keep],mag_list[keep],'k.',markersize=0.75)
	axi.plot(cent_mark, mean_val,'ro',markersize=10)
	axi.set_ylim((0,10))
	axi.set_xlabel('distance to cell surface')
	axi.set_ylabel(r'displacement magnitude $\mu m$')

# --> 3D plot of the cell, configuration number influences title and color 
def plot_cell_3D(ax,cell_num,cell_mesh, cell_center, cell_vol, X_DIM, Y_DIM, Z_DIM):
	if cell_num == 1:
		col = (0.75,0.75,0.75)
	elif cell_num == 2:
		col = (0,0,0)
	verts = cell_mesh; cent = cell_center; vol = cell_vol
	
	ax.set_aspect('auto')
	ax.plot(verts[:,0],verts[:,1],verts[:,2],'.',color=col)
	
	ax.set_xlim((-1,X_DIM))
	ax.set_ylim((-1,Y_DIM))
	ax.set_zlim((-1,Z_DIM))
		
	if cell_num == 1:
		ax.set_title('cell config 1, %.1f $\mu m^3$'%(vol))
	elif cell_num == 2:
		ax.set_title('cell config 2, %.1f $\mu m^3$'%(vol))

# --> plot of scores (type 1 is similarity to neighbors, type 2 is direction relative to cell)
def plot_scores_subplot(data,title,axi,color_type):
	num_pts = 250
	X_plot = np.linspace(np.min(data),np.max(data),num_pts).reshape(-1,1)
	X =  data.reshape(-1,1)
	kde = KernelDensity(kernel='gaussian', bandwidth=0.1).fit(X)
	log_dens = kde.score_samples(X_plot)
	axi.set_xlabel('score')
	axi.set_ylabel('probability density function')
	axi.set_title(title)
	ci_max = np.max(data); ci_min = np.min(data) 
	axi.plot(X_plot[:,0],np.exp(log_dens),'k-',linewidth=0.5)
	for kk in range(0,num_pts):
		ci = X_plot[kk,0]
		if color_type == 1: #--> all positive numbers, blue is 0, red is high 
			col = ((ci-ci_min)/(ci_max-ci_min), 0, 1.0 - (ci-ci_min)/(ci_max-ci_min))
		elif color_type == 2:
			if ci < 0.0:
				col = (np.abs(ci),0,0.5*np.abs(ci))
			else:
				col = (0, np.abs(ci), np.abs(ci))
		axi.plot(X_plot[kk, 0], np.exp(log_dens[kk]),'.',color=col)
	return

# --> helper function plots slice of cell
def plot_cell(cent,project_1,project_2,project_out,col,axi):
	buffer = 0.5
	buffer_up = cent + buffer
	buffer_low = cent - buffer
	plot_1 = []
	plot_2 = [] 
	num_pts = project_1.shape[0]
	for kk in range(0,num_pts):
		if project_out[kk] < buffer_up and project_out[kk] > buffer_low:
			plot_1.append(project_1[kk])
			plot_2.append(project_2[kk])
	axi.plot(plot_1,plot_2,'.',color=col)
	return

# --> helper function plots slice of vectors 
def plot_vectors(color_type, color_info, project_1, project_2, project_1d, project_2d, cent, project_out, axi):
	ci_min = np.min(color_info); ci_max = np.max(color_info)
	num_pts = project_1.shape[0]
	for kk in range(0,num_pts):
		# --> the vectors themselves 
		scale = 1 
		proj1_a = project_1[kk]; proj1_d = project_1d[kk]*scale
		proj2_a = project_2[kk]; proj2_d = project_2d[kk]*scale
		pout = project_out[kk]; 
		buffer = 10 
		# --> color of the vectors 
		ci = color_info[kk] 
		
		if pout > cent - buffer and pout < cent + buffer:
			# --> colortype 
			if color_type == 1: #--> all positive numbers, blue is 0, red is high 
				col = ((ci-ci_min)/(ci_max-ci_min), 0, 1.0 - (ci-ci_min)/(ci_max-ci_min))
			elif color_type == 2:
				if ci < 0.0:
					col = (np.abs(ci),0,0.5*np.abs(ci))
				else:
					col = (0, np.abs(ci), np.abs(ci))
			# --> plot the vectors 
			axi.arrow(proj1_a,proj2_a,proj1_d,proj2_d,color = col,linewidth=1.0,head_width=1.5)
	return 

# --> plot a slice plot, each has beads and a cell 
def plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM):
	num_beads = X.shape[0]
	XYZ = np.zeros((num_beads,3)); XYZ[:,0] = X; XYZ[:,1] = Y; XYZ[:,2] = Z 
	UVW = np.zeros((num_beads,3)); UVW[:,0] = U; UVW[:,1] = V; UVW[:,2] = W
	
	cell_center_avg = 0.5*cell_center_1 + 0.5*cell_center_2
	
	if plane_type == 1: #XZ-plane
		idx_1 = 0
		idx_2 = 2
		idx_out = 1
	elif plane_type == 2: #YZ-plane
		idx_1 = 1
		idx_2 = 2
		idx_out = 0
	elif plane_type == 3: #XY-plane 
		idx_1 = 0
		idx_2 = 1
		idx_out = 2 
		
	cent = cell_center_avg[idx_out]
	project_1_cell_A = cell_mesh_1[:,idx_1]
	project_2_cell_A = cell_mesh_1[:,idx_2]
	project_out_cell_A = cell_mesh_1[:,idx_out]
	cell_color_A = (0.75,0.75,0.75)
	
	project_1_cell_B = cell_mesh_2[:,idx_1]
	project_2_cell_B = cell_mesh_2[:,idx_2]
	project_out_cell_B = cell_mesh_2[:,idx_out]
	cell_color_B = (0.0,0.0,0.0)
	
	project_1_bead = XYZ[:,idx_1]
	project_2_bead = XYZ[:,idx_2]
	project_1d_bead = UVW[:,idx_1]
	project_2d_bead = UVW[:,idx_2]
	project_out_bead = XYZ[:,idx_out]
	
	# call cell plot for cell 1 
	plot_cell(cent,project_1_cell_A,project_2_cell_A,project_out_cell_A,cell_color_A,axi)
	# call cell plot for cell 2 
	plot_cell(cent,project_1_cell_B,project_2_cell_B,project_out_cell_B,cell_color_B,axi)
	# call vector plot 
	plot_vectors(color_type, color_info, project_1_bead, project_2_bead, project_1d_bead, project_2d_bead, cent, project_out_bead, axi)
	
	center = cell_center_avg
	
	if plane_type == 1: #XZ-plane
		axi.plot([-1,X_DIM],[center[2],center[2]],'k:',linewidth=1.0)
		axi.plot([center[0],center[0]],[-1,Z_DIM],'k:',linewidth=1.0)
		axi.set_xlim((-1,X_DIM)) 
		axi.set_ylim((-1,Z_DIM))
		axi.set_xlabel(r'x-position $\mu m$')
		axi.set_ylabel(r'z-position $\mu m$')
	elif plane_type == 2: #YZ-plane
		axi.plot([-1,Y_DIM],[center[2],center[2]],'k:',linewidth=1.0)
		axi.plot([center[1],center[1]],[-1,Z_DIM],'k:',linewidth=1.0)
		axi.set_xlim((-1,Y_DIM)) 
		axi.set_ylim((-1,Z_DIM))
		axi.set_xlabel(r'y-position $\mu m$')
		axi.set_ylabel(r'z-position $\mu m$')
	elif plane_type == 3: #XY-plane
		axi.plot([-1,X_DIM],[center[1],center[1]],'k:',linewidth=1.0)
		axi.plot([center[0],center[0]],[-1,Y_DIM],'k:',linewidth=1.0)
		axi.set_xlim((-1,X_DIM)) 
		axi.set_ylim((-1,Y_DIM))
		axi.set_xlabel(r'x-position $\mu m$')
		axi.set_ylabel(r'y-position $\mu m$')
	return 

# --> plot a cell-vector row 
def plot_cell_vector_slice_row(ax_list,color_type,color_info,X,Y,Z,U,V,W,cell_center_1,cell_mesh_1,cell_center_2,cell_mesh_2,X_DIM,Y_DIM,Z_DIM):
	axi = ax_list[0] 
	plane_type = 1 
	plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = ax_list[1]
	plane_type = 2
	plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = ax_list[2] 
	plane_type = 3
	plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	return 

# --> plot cells
def plot_only_cells(cell_mesh_1,cell_center_1,cell_vol_1,cell_mesh_2,cell_center_2,cell_vol_2,X_DIM,Y_DIM,Z_DIM,folder,figtype_list):
	fig = plt.figure()
	plt.style.use(stylepath)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig.set_figheight(5)
	fig.set_figwidth(10)
	axi = fig.add_subplot(1, 2, 1, projection='3d')
	plot_cell_3D(axi,1,cell_mesh_1, cell_center_1, cell_vol_1, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(1, 2, 2, projection='3d')
	plot_cell_3D(axi,2,cell_mesh_2, cell_center_2, cell_vol_2, X_DIM, Y_DIM, Z_DIM)
	plt.tight_layout()
	for end in figtype_list:
		fname = folder + '/Cell_plots_3D' + end
		plt.savefig(fname)
	return  

# --> plot scores
def plot_only_scores(neigh_score,dir_score,folder,figtype_list):
	fig = plt.figure()
	plt.style.use(stylepath)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig.set_figheight(5)
	fig.set_figwidth(10)
	axi = fig.add_subplot(1,2,1)
	data = np.asarray(neigh_score)
	color_type = 1 
	title = 'neighbor distance score'
	plot_scores_subplot(data,title,axi,color_type)
	axi = fig.add_subplot(1,2,2)
	data = np.asarray(dir_score)
	color_type = 2
	title = r'$n_{cell} \cdot n_{vector}$'
	plot_scores_subplot(data,title,axi,color_type)
	plt.tight_layout()
	for end in figtype_list:
		fname = folder + '/Score_plots' + end
		plt.savefig(fname)
	return

# --> plot slice
def plot_only_slice(dir_score,X,Y,Z,U,V,W,cell_center_1,cell_mesh_1,cell_center_2,cell_mesh_2,X_DIM,Y_DIM,Z_DIM,folder,figtype_list):
	fig = plt.figure()
	plt.style.use(stylepath)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig.set_figheight(5)
	fig.set_figwidth(15)
	color_type = 2
	color_info = dir_score
	axi = fig.add_subplot(1,3,1)
	plane_type = 1
	plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(1,3,2)
	plane_type = 2
	plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(1,3,3)
	plane_type = 3
	plot_cell_vector_slice(color_type, color_info, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	plt.tight_layout()
	for end in figtype_list:
		fname = folder + '/Bead_disp_slice' + end
		plt.savefig(fname)
	return

# --> plot distance
def plot_only_distance(cell_mesh,dist_from_edge,dist_from_cell,mag_list,folder,figtype_list):
	fig = plt.figure()
	plt.style.use(stylepath)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig.set_figheight(5)
	fig.set_figwidth(5)
	axi = fig.add_subplot(1,1,1)
	plot_surface_disp(axi,cell_mesh,dist_from_edge,dist_from_cell, mag_list)
	plt.tight_layout()
	for end in figtype_list:
		fname = folder + '/Disp_wrt_dist' + end
		plt.savefig(fname)
	return 

# --> plot all
def plot_all(folder, root_directory, file_prefix_1,file_prefix_2,dir_score,neigh_score,dist_from_edge,dist_from_cell,mag_list,\
		X,Y,Z,U,V,W,cell_center_1,cell_mesh_1,cell_vol_1,cell_center_2,cell_mesh_2,cell_vol_2,X_DIM,Y_DIM,Z_DIM,figtype_list):
	fig = plt.figure()
	plt.style.use(stylepath)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig.set_figheight(10)
	fig.set_figwidth(20)
	axi = fig.add_subplot(2,4,1)
	data = np.asarray(dir_score)
	color_type = 2
	title = r'$n_{cell} \cdot n_{vector}$'
	plot_scores_subplot(data,title,axi,color_type)
	axi = fig.add_subplot(2,4,2)
	plane_type = 1
	plot_cell_vector_slice(color_type, dir_score, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(2,4,3)
	plane_type = 2
	plot_cell_vector_slice(color_type, dir_score, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(2,4,4)
	plane_type = 3
	plot_cell_vector_slice(color_type, dir_score, X, Y, Z, U, V, W, cell_center_1,\
		cell_mesh_1, cell_center_2, cell_mesh_2, plane_type, axi, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(2,4,5)
	data = np.asarray(neigh_score)
	color_type = 1 
	title = 'neighbor distance score'
	plot_scores_subplot(data,title,axi,color_type)
	axi = fig.add_subplot(2,4,6)
	plot_surface_disp(axi,cell_mesh_1,dist_from_edge,dist_from_cell, mag_list)
	axi = fig.add_subplot(2,4,7, projection='3d')
	plot_cell_3D(axi,1,cell_mesh_1, cell_center_1, cell_vol_1, X_DIM, Y_DIM, Z_DIM)
	axi = fig.add_subplot(2,4,8, projection='3d')
	plot_cell_3D(axi,2,cell_mesh_2, cell_center_2, cell_vol_2, X_DIM, Y_DIM, Z_DIM)
	
	plt.tight_layout()
	for end in figtype_list:
		fname = folder + '/Summary_plot' + end
		plt.savefig(fname)
	for end in figtype_list:
		fname = root_directory + '/Post_proc_summary' + '/' + 'Summary_' + file_prefix_1 + '_to_' + file_prefix_2 + end
		plt.savefig(fname) 
	return 

# call individual plots, plus call multiple subplots
def call_plot_main(plot_type,file_prefix_1,file_prefix_2,num_feat,X_DIM,Y_DIM,Z_DIM,figtype_list,use_corrected_cell,root_directory):
	folder = root_directory + '/Track_' + file_prefix_1 + '_to_' + file_prefix_2 
	if use_corrected_cell:
		cell_mesh_2 = np.loadtxt(folder + '/cell_mesh_2_corrected.txt')
	X, Y, Z, U, V, W = import_bead_disps(folder)
	cell_mesh_1, cell_normal_1, cell_center_1, cell_vol_1, cell_mesh_2, cell_normal_2, cell_center_2, cell_vol_2 = import_cell_info(file_prefix_1,file_prefix_2,root_directory)
	
	neigh_score = color_point_neighbor_similarity(X, Y, Z, U, V, W, num_feat)
	dir_score, dist_from_cell, mag_list = color_point_direction(X, Y, Z, U, V, W, cell_mesh_1, cell_normal_1)
	dist_from_edge = compute_dist_from_edge(X, Y, Z, X_DIM, Y_DIM, Z_DIM)
	
	#type 6 will create all plots
	# --> arrange data 
	if plot_type == 1 or plot_type == 6: # big plot with everything, saves it in two directories
		plot_all(folder, root_directory, file_prefix_1,file_prefix_2,dir_score,neigh_score,dist_from_edge,dist_from_cell,mag_list,\
			X,Y,Z,U,V,W,cell_center_1,cell_mesh_1,cell_vol_1,cell_center_2,cell_mesh_2,cell_vol_2,X_DIM,Y_DIM,Z_DIM,figtype_list)
	if plot_type == 2 or plot_type == 6: # plots cells in both configurations 
		plot_only_cells(cell_mesh_1,cell_center_1,cell_vol_1,cell_mesh_2,cell_center_2,cell_vol_2,X_DIM,Y_DIM,Z_DIM,folder,figtype_list)
	if plot_type == 3 or plot_type == 6: # plots scores only 
		plot_only_scores(neigh_score,dir_score,folder,figtype_list)
	if plot_type == 4 or plot_type == 6: # plots slice only 
		plot_only_slice(dir_score,X,Y,Z,U,V,W,cell_center_1,cell_mesh_1,cell_center_2,cell_mesh_2,X_DIM,Y_DIM,Z_DIM,folder,figtype_list)
	if plot_type == 5 or plot_type == 6: # plots magnitude wrt distance from surface
		plot_only_distance(cell_mesh_1,dist_from_edge,dist_from_cell,mag_list,folder,figtype_list)
	return

##########################################################################################
# displacement interpolation -- use GPR 
##########################################################################################
# --> create GP model 
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

# --> create GP models 
def create_GP_model(file_prefix_1,file_prefix_2,root_directory):
	folder = root_directory + '/Track_' + file_prefix_1 + '_to_' + file_prefix_2
	X, Y, Z, U, V, W = import_bead_disps(folder)
	gp_U, scaler = create_gp_model(X,Y,Z,U)
	gp_V, scaler = create_gp_model(X,Y,Z,V)
	gp_W, scaler = create_gp_model(X,Y,Z,W)
	pickle.dump(gp_U, open(folder + '/gp_U.sav','wb'))
	pickle.dump(gp_V, open(folder + '/gp_V.sav','wb'))
	pickle.dump(gp_W, open(folder + '/gp_W.sav','wb'))
	pickle.dump(scaler,open(folder + '/scaler.sav','wb'))
	return 

# --> interpolate GP model 
def interpolate_gp_model(plane_case, center, gp, scaler, X_DIM, Y_DIM, Z_DIM ):
	x_min = -1; x_max = X_DIM
	y_min = -1; y_max = Y_DIM
	z_min = -1; z_max = Z_DIM
	grid_pts = 100 
	# --> construct artificial grid for plotting 
	if plane_case == 1: #x plane
		x = center[1]
		y = np.linspace(y_min,y_max,grid_pts)
		z = np.linspace(z_min,z_max,grid_pts)
		Y, Z = np.meshgrid(y,z)
		X = x * np.ones((grid_pts,grid_pts))
		RES = np.zeros((grid_pts,grid_pts))
	elif plane_case == 2: #y plane
		x = np.linspace(x_min,x_max,grid_pts)
		y = center[0]
		z = np.linspace(z_min,z_max,grid_pts)
		X, Z = np.meshgrid(x, z)
		Y = y * np.ones((grid_pts,grid_pts))
		RES = np.zeros((grid_pts,grid_pts))
	elif plane_case == 3: #z plane
		x = np.linspace(x_min,x_max,grid_pts)
		y = np.linspace(y_min,y_max,grid_pts)
		z = center[2]
		X, Y = np.meshgrid(x, y)
		Z = z * np.ones((grid_pts,grid_pts))
		RES = np.zeros((grid_pts,grid_pts))
	
	# --> fit model grid 
	for j in range(0,grid_pts):
		input = []
		for k in range(0,grid_pts):
			li = [X[j,k],Y[j,k],Z[j,k]]
			input.append(li)
		input = np.asarray(input)
		input = scaler.transform(input)
		pred = gp.predict(input)
		RES[j,:] = pred[:]
	
	if plane_case == 1:
		return Y, Z, RES
	elif plane_case == 2:
		return X, Z, RES
	elif plane_case == 3:
		return X, Y, RES

# --> create a single GP plot 
def plot_gp_model_single_plot(axi,is_mag,data_1,data_2,result,title):
	# --> plot interpolated field
	vmin = -5; vmax = 5
	if is_mag:
		vmin = 0; vmax = 10 
	CS1 = axi.pcolor(data_1, data_2, result, cmap=plt.cm.coolwarm,vmin=vmin,vmax=vmax)
	cbar = plt.colorbar(CS1, ax=axi)
	cbar.set_label(title,labelpad=-95,y=1.13,rotation=0)
	return

# --> plot GPR model, one row 
def plot_gp_model_one_row(ax_list,is_mag,X_DIM,Y_DIM,Z_DIM,title,center,gp_model,scaler,cell_mesh_1,cell_mesh_2):	
	axi = ax_list[0]
	plane_case = 2
	if is_mag == False:
		data_1, data_2, result = interpolate_gp_model(plane_case, center, gp_model, scaler, X_DIM, Y_DIM, Z_DIM)
	else:
		data_1, data_2, result_0 = interpolate_gp_model(plane_case, center, gp_model[0], scaler, X_DIM, Y_DIM, Z_DIM)
		data_1, data_2, result_1 = interpolate_gp_model(plane_case, center, gp_model[1], scaler, X_DIM, Y_DIM, Z_DIM)
		data_1, data_2, result_2 = interpolate_gp_model(plane_case, center, gp_model[2], scaler, X_DIM, Y_DIM, Z_DIM)
		result = (result_0**2.0 + result_1**2.0 + result_2**2.0)**(1.0/2.0)
	plot_gp_model_single_plot(axi,is_mag,data_1,data_2,result,title)
	idx0 = 0; idx1 = 2; idx2 = 1
	plot_cell(center[idx2],cell_mesh_1[:,idx0],cell_mesh_1[:,idx1],cell_mesh_1[:,idx2],(0.75,0.75,0.75),axi)
	plot_cell(center[idx2],cell_mesh_2[:,idx0],cell_mesh_2[:,idx1],cell_mesh_2[:,idx2],(0,0,0),axi)
	axi.plot([-1,X_DIM],[center[2],center[2]],'k:',linewidth=1.0)
	axi.plot([center[0],center[0]],[-1,Z_DIM],'k:',linewidth=1.0)
	axi.set_xlabel(r'x-position $\mu m$')
	axi.set_ylabel(r'z-position $\mu m$')
	axi.set_xlim((-1,X_DIM)) 
	axi.set_ylim((-1,Z_DIM))
	
	axi = ax_list[1]
	place_case = 1 
	if is_mag == False:
		data_1, data_2, result = interpolate_gp_model(plane_case, center, gp_model, scaler, X_DIM, Y_DIM, Z_DIM)
	else:
		data_1, data_2, result_0 = interpolate_gp_model(plane_case, center, gp_model[0], scaler, X_DIM, Y_DIM, Z_DIM)
		data_1, data_2, result_1 = interpolate_gp_model(plane_case, center, gp_model[1], scaler, X_DIM, Y_DIM, Z_DIM)
		data_1, data_2, result_2 = interpolate_gp_model(plane_case, center, gp_model[2], scaler, X_DIM, Y_DIM, Z_DIM)
		result = (result_0**2.0 + result_1**2.0 + result_2**2.0)**(1.0/2.0)
	plot_gp_model_single_plot(axi,is_mag,data_1,data_2,result,title)
	idx0 = 1; idx1 = 2; idx2 = 0 
	plot_cell(center[idx2],cell_mesh_1[:,idx0],cell_mesh_1[:,idx1],cell_mesh_1[:,idx2],(0.75,0.75,0.75),axi)
	plot_cell(center[idx2],cell_mesh_2[:,idx0],cell_mesh_2[:,idx1],cell_mesh_2[:,idx2],(0,0,0),axi)
	axi.plot([-1,Y_DIM],[center[2],center[2]],'k:',linewidth=1.0)
	axi.plot([center[1],center[1]],[-1,Z_DIM],'k:',linewidth=1.0)
	axi.set_xlabel(r'y-position $\mu m$')
	axi.set_ylabel(r'z-position $\mu m$')
	axi.set_xlim((-1,Y_DIM)) 
	axi.set_ylim((-1,Z_DIM))
	
	axi = ax_list[2]
	plane_case = 3
	if is_mag == False:
		data_1, data_2, result = interpolate_gp_model(plane_case, center, gp_model, scaler, X_DIM, Y_DIM, Z_DIM)
	else:
		data_1, data_2, result_0 = interpolate_gp_model(plane_case, center, gp_model[0], scaler, X_DIM, Y_DIM, Z_DIM)
		data_1, data_2, result_1 = interpolate_gp_model(plane_case, center, gp_model[1], scaler, X_DIM, Y_DIM, Z_DIM)
		data_1, data_2, result_2 = interpolate_gp_model(plane_case, center, gp_model[2], scaler, X_DIM, Y_DIM, Z_DIM)
		result = (result_0**2.0 + result_1**2.0 + result_2**2.0)**(1.0/2.0)
	plot_gp_model_single_plot(axi,is_mag,data_1,data_2,result,title)
	idx0 = 0; idx1 = 1; idx2 = 2
	plot_cell(center[idx2],cell_mesh_1[:,idx0],cell_mesh_1[:,idx1],cell_mesh_1[:,idx2],(0.75,0.75,0.75),axi)
	plot_cell(center[idx2],cell_mesh_2[:,idx0],cell_mesh_2[:,idx1],cell_mesh_2[:,idx2],(0,0,0),axi)
	axi.plot([-1,X_DIM],[center[1],center[1]],'k:',linewidth=1.0)
	axi.plot([center[0],center[0]],[-1,Z_DIM],'k:',linewidth=1.0)
	axi.set_xlabel(r'x-position $\mu m$')
	axi.set_ylabel(r'y-position $\mu m$')	
	axi.set_xlim((-1,X_DIM)) 
	axi.set_ylim((-1,Z_DIM))
	

# --> plot GPR model 
def plot_gp_model(file_prefix_1,file_prefix_2,X_DIM,Y_DIM,Z_DIM,figtype_list,use_corrected_cell, root_directory):
	cell_mesh_1, cell_normal_1, cell_center_1, cell_vol_1, cell_mesh_2, cell_normal_2, cell_center_2, cell_vol_2 = import_cell_info(file_prefix_1,file_prefix_2)
	
	center = 0.5*cell_center_1 + 0.5*cell_center_2
	folder = root_directory + '/Track_' + file_prefix_1 + '_to_' + file_prefix_2
	
	if use_corrected_cell:
		cell_mesh_2 = np.loadtxt(folder + '/cell_mesh_2_corrected.txt')
		
	gp_U = pickle.load(open(folder + '/gp_U.sav', 'rb'))
	gp_V = pickle.load(open(folder + '/gp_V.sav', 'rb'))
	gp_W = pickle.load(open(folder + '/gp_W.sav', 'rb'))
	scaler = pickle.load(open(folder + '/scaler.sav','rb'))
	
	fig = plt.figure()
	plt.style.use(stylepath)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig.set_figheight(20)
	fig.set_figwidth(15)
	
	ax1 = fig.add_subplot(4, 3, 1); ax2 = fig.add_subplot(4, 3, 2);ax3 = fig.add_subplot(4, 3, 3) 
	ax_list1 = [ax1,ax2,ax3]
	title = r'x-displacement $\mu m$'
	is_mag = False
	plot_gp_model_one_row(ax_list1,is_mag,X_DIM,Y_DIM,Z_DIM,title,center,gp_U,scaler,cell_mesh_1,cell_mesh_2)
	
	ax4 = fig.add_subplot(4, 3, 4); ax5 = fig.add_subplot(4, 3, 5); ax6 = fig.add_subplot(4, 3, 6) 
	ax_list2 = [ax4,ax5,ax6]
	title = r'y-displacement $\mu m$'
	is_mag = False
	plot_gp_model_one_row(ax_list2,is_mag,X_DIM,Y_DIM,Z_DIM,title,center,gp_V,scaler,cell_mesh_1,cell_mesh_2)
	
	ax7 = fig.add_subplot(4, 3, 7); ax8 = fig.add_subplot(4, 3, 8); ax9 = fig.add_subplot(4, 3, 9) 
	ax_list3 = [ax7,ax8,ax9]
	title = r'z-displacement $\mu m$'
	is_mag = False
	plot_gp_model_one_row(ax_list3,is_mag,X_DIM,Y_DIM,Z_DIM,title,center,gp_W,scaler,cell_mesh_1,cell_mesh_2)
	
	ax10 = fig.add_subplot(4, 3, 10); ax11 = fig.add_subplot(4, 3, 11); ax12 = fig.add_subplot(4, 3, 12) 
	ax_list4 = [ax10,ax11,ax12]
	title = r'mag-displacement $\mu m$'
	is_mag = True
	plot_gp_model_one_row(ax_list4,is_mag,X_DIM,Y_DIM,Z_DIM,title,center,[gp_U, gp_V, gp_W],scaler,cell_mesh_1,cell_mesh_2)
	
	plt.tight_layout()
	for end in figtype_list:
		fname = folder + '/Interpolate_plot' + end
		plt.savefig(fname)