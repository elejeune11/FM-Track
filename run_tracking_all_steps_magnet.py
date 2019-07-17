import input_info_magnet as info
import post_process as post_process 
import pre_process as pre_process
import tracking as track 
import numpy as np
import os 
##########################################################################################
# indicate which steps to run
##########################################################################################
run_pre_process = True
run_tracking = True
run_post_process = True
##########################################################################################
# get set up info from input files 
##########################################################################################
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
num_tracking_pairs = len(tracking_pairs)
if run_tracking:
	for kk in range(0,num_tracking_pairs):
		closest_no_conflict = track.track_main_call(track_type,tracking_pairs[kk][0],tracking_pairs[kk][1], num_feat,num_nearest, buffer_cell)
##########################################################################################
# post process 
##########################################################################################
num_tracking_pairs = len(tracking_pairs)
if run_post_process:
	for kk in range(0,num_tracking_pairs):
		post_process.call_plot_main(plot_type,tracking_pairs[kk][0], tracking_pairs[kk][1],num_feat,\
			X_DIM,Y_DIM,Z_DIM,figtype_list, use_corrected_cell)
		
	if run_GP:
		for kk in range(0,num_tracking_pairs):
			post_process.create_GP_model(tracking_pairs[kk][0], tracking_pairs[kk][1])
			post_process.plot_gp_model(tracking_pairs[kk][0], tracking_pairs[kk][1],\
				X_DIM,Y_DIM,Z_DIM,figtype_list)