from . import post_process
from . import pre_process
from . import tracking as track
import numpy as np
import os 

##########################################################################################
# get set up info from input files 
##########################################################################################

def run_tracking_all_steps(run_pre_process, run_tracking, run_post_process, info):
	filenames_beads, dirnames_beads = info.get_filenames_beads()
	cell_channel, bead_channel = info.get_color_channels()
	X_DIM, Y_DIM, Z_DIM = info.get_FOV_dims(3)
	cell_threshold = info.get_cell_thresh()
	filenames_cell, dirnames_cell  = info.get_filenames_cell()
	savefnames_cell, savefnames_beads  = info.get_savenames()
	tracking_pairs = info.get_tracking_pairs()
	num_feat, num_nearest, buffer_cell, track_type  = info.get_tracking_params()
	figtype_list, plot_type, run_GP, use_corrected_cell = info.get_postproc_info()
	##########################################################################################
	# pre-process
	##########################################################################################
	print('Running Pre-Processing Step')
	num_imgs = len(filenames_beads)
	if run_pre_process:
		for kk in range(0,num_imgs):
			pre_process.get_bead_centers(dirnames_beads[kk],filenames_beads[kk],\
				savefnames_beads[kk],bead_channel, X_DIM, Y_DIM, Z_DIM)
			pre_process.get_cell_surface(dirnames_cell[kk],filenames_cell[kk],\
				savefnames_cell[kk],cell_channel, X_DIM, Y_DIM, Z_DIM, cell_threshold)
	##########################################################################################
	# run tracking 
	##########################################################################################
	print('Tracking')
	num_tracking_pairs = len(tracking_pairs)
	if run_tracking:
		for kk in range(0,num_tracking_pairs):
			closest_no_conflict = track.track_main_call(track_type,tracking_pairs[kk][0],tracking_pairs[kk][1], num_feat,num_nearest, buffer_cell, info.root_directory)
	##########################################################################################
	# post process 
	##########################################################################################
	print('Running Post-Processing Step')
	num_tracking_pairs = len(tracking_pairs)
	if run_post_process:
		for kk in range(0,num_tracking_pairs):
			post_process.call_plot_main(plot_type,tracking_pairs[kk][0], tracking_pairs[kk][1],num_feat,\
				X_DIM,Y_DIM,Z_DIM,figtype_list, use_corrected_cell, info.root_directory)
			
		if run_GP:
			for kk in range(0,num_tracking_pairs):
				post_process.create_GP_model(tracking_pairs[kk][0], tracking_pairs[kk][1], info.root_directory)
				post_process.plot_gp_model(tracking_pairs[kk][0], tracking_pairs[kk][1],\
					X_DIM,Y_DIM,Z_DIM,figtype_list, info.root_directory)