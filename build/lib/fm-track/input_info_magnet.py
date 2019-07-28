##########################################################################################
# import packages
##########################################################################################
import os 
##########################################################################################
# functions to get path names and constants that have to be hard coded in 
##########################################################################################
def get_filenames_cell():
	filenames_cell = [ \
		'data/7_16_Magnet/Magnet/Magnetic Bead/Magnet%s.tif',\
		'data/7_16_Magnet/No Magnet/Magnetic Bead/No Magnet%s.tif',\
		'data/7_16_Magnet/No Magnet 2/Magnetic Bead/No Magnet 2%s.tif']
	dirnames_cell = [ \
		'data/7_16_Magnet/Magnet/Magnetic Bead/',\
		'data/7_16_Magnet/No Magnet/Magnetic Bead/',\
		'data/7_16_Magnet/No Magnet 2/Magnetic Bead/']
	return filenames_cell, dirnames_cell 
	
def get_filenames_beads():
	filenames_beads = [ \
		'data/7_16_Magnet/Magnet/Tracking Beads/Magnet%s.tif',\
		'data/7_16_Magnet/No Magnet/Tracking Beads/No Magnet%s.tif',\
		'data/7_16_Magnet/No Magnet 2/Tracking Beads/No Magnet 2%s.tif']
	dirnames_beads = [ \
		'data/7_16_Magnet/Magnet/Tracking Beads/',\
		'data/7_16_Magnet/No Magnet/Tracking Beads/',\
		'data/7_16_Magnet/No Magnet 2/Tracking Beads/']
	return filenames_beads, dirnames_beads
	
def get_savenames():
	out_folder_cell = 'Gel_cell_coords'
	if not os.path.exists(out_folder_cell):
		os.makedirs(out_folder_cell)
	out_folder_beads = 'Gel_bead_center_coords'
	if not os.path.exists(out_folder_beads):
		os.makedirs(out_folder_beads)
	savefnames = [ \
		'data/7_16_Magnet',\
		'data/7_16_NoMagnet',\
		'data/7_16_NoMagnet2']
	savefnames_cell = [] 
	savefnames_beads = [] 
	for kk in range(0,len(savefnames)):
		savefnames_cell.append(out_folder_cell + '/' + savefnames[kk] + '_cell_')
		savefnames_beads.append(out_folder_beads + '/' + savefnames[kk] + '_beads.txt') 
	return savefnames_cell, savefnames_beads 

def get_tracking_pairs():
	out_folder = 'Post_proc_summary'
	if not os.path.exists(out_folder):
		os.makedirs(out_folder)
	tracking_pairs = [ \
		['data/7_16_Magnet','../7_16_NoMagnet'],\
		['data/7_16_NoMagnet','../7_16_NoMagnet2'],\
		['data/7_16_Magnet','../7_16_NoMagnet2']]
	return tracking_pairs

def get_color_channels(): #for indexing 3D array with cell image and bead image
	cell_channel = 0 #CellBrite Red
	bead_channel = 1 #Green fluorescent beads 
	return cell_channel, bead_channel

def get_FOV_dims(type):
	if type == 1: # obj_63x
		X_DIM = 149.95
		Y_DIM = 149.95
		Z_DIM = 140.0
	elif type == 2: # obj_???
		X_DIM = 141.70
		Y_DIM = 141.70
		Z_DIM = 120.0 
	elif type == 3:
		X_DIM = 149.95
		Y_DIM = 149.95
		Z_DIM = 120.0
	return X_DIM, Y_DIM, Z_DIM 
	
def get_cell_thresh():
	cell_thresh = 1.0
	return cell_thresh
	
def get_tracking_params():
	num_feat = 5
	num_nearest = 15
	buffer_cell = 0
	track_type = 2 # type 1 will NOT perform translation correction, type 2 will
	return num_feat, num_nearest, buffer_cell, track_type 

def get_translation_correction_params():
	buffer_cell = 30
	return buffer_cell 

def get_postproc_info():
	figtype_list = ['.png'] 
	plot_type = 6.0
	run_GP = False
	use_corrected_cell = True
	return figtype_list, plot_type, run_GP, use_corrected_cell

def get_list_of_packages():
	#list of all packages needed to run this code 
	package_list = [] 
	return package_list 