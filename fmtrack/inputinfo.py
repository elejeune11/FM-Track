##########################################################################################
# import packages
##########################################################################################
import os 
import numpy as np
from . import pre_process
from . import post_process
from . import tracking as track

##########################################################################################
# class to get path names and constants that can be modified
##########################################################################################

class InputInfo:

	def __init__(self,root_directory=None):

		self.root_directory = root_directory
		self.assign_defaults()

	def assign_defaults(self):

		# For get_color_channels()
		self.cell_channel = 0 #CellBrite Red
		self.bead_channel = 1 #Green fluorescent beads 

		# For get_cell_thresh()
		self.cell_thresh = 1.0

		# For get_tracking_params()
		self.num_feat = 5
		self.num_nearest = 15
		self.buffer_cell_tracking = 0
		self.track_type = 2 # type 1 will NOT perform translation correction, type 2 will

		# For get_translation_correction_params()
		self.buffer_cell_translation = 30

		# For get_postproc_info()
		self.figtype_list = ['.png'] 
		self.plot_type = 6.0
		self.run_GP = False
		self.use_corrected_cell = True
		self.should_plot = True
		self.print_progress = True

		self.set_out_folder_cell('Gel_cell_coords')
		self.set_out_folder_beads('Gel_bead_center_coords')
		self.set_out_folder('Post_proc_summary')

	def run_tracking_all_steps(self,run_pre_process, run_tracking, run_post_process):
		if run_pre_process:
			self.run_pre_process()
		if run_tracking:
			self.run_tracking()
		if run_post_process:
			self.run_post_process()

	def run_pre_process(self):
		if self.print_progress:
			print('Running Pre-Processing Step')
		filenames_beads, dirnames_beads = self.get_filenames_beads()
		filenames_cell, dirnames_cell  = self.get_filenames_cell()
		savefnames_cell, savefnames_beads  = self.get_savenames()
		cell_channel, bead_channel = self.get_color_channels()
		cell_threshold = self.get_cell_thresh()
		X_DIM, Y_DIM, Z_DIM = self.get_FOV_dims()
		num_imgs = len(filenames_beads)
		self.bead_array = []
		self.mesh_array = []
		for kk in range(0,num_imgs):
			self.bead_array.append(pre_process.get_bead_centers(dirnames_beads[kk],filenames_beads[kk],\
				savefnames_beads[kk],bead_channel, X_DIM, Y_DIM, Z_DIM))
			self.mesh_array.append(pre_process.get_cell_surface(dirnames_cell[kk],filenames_cell[kk],\
				savefnames_cell[kk],cell_channel, X_DIM, Y_DIM, Z_DIM, cell_threshold))

	def run_tracking(self):
		if self.print_progress:
			print('Tracking:')
		num_feat, num_nearest, buffer_cell, track_type  = self.get_tracking_params()
		tracking_pairs = self.get_tracking_pairs()
		num_tracking_pairs = len(tracking_pairs)
		for kk in range(0,num_tracking_pairs):
			_ = track.track_main_call(track_type,tracking_pairs[kk][0],tracking_pairs[kk][1], num_feat,num_nearest, buffer_cell, self)

	def run_post_process(self):
		if self.print_progress:
			print('Running Post-Processing Step')
		X_DIM, Y_DIM, Z_DIM = self.get_FOV_dims()
		figtype_list, plot_type, run_GP, use_corrected_cell = self.get_postproc_info()
		num_feat, num_nearest, buffer_cell, track_type  = self.get_tracking_params()
		tracking_pairs = self.get_tracking_pairs()
		num_tracking_pairs = len(tracking_pairs)
		for kk in range(0,num_tracking_pairs):
			post_process.call_plot_main(plot_type,tracking_pairs[kk][0], tracking_pairs[kk][1],num_feat,\
				X_DIM,Y_DIM,Z_DIM,figtype_list, use_corrected_cell, self.root_directory, self.should_plot)
			
		if run_GP:
			for kk in range(0,num_tracking_pairs):
				post_process.create_GP_model(tracking_pairs[kk][0], tracking_pairs[kk][1], self.root_directory)
				post_process.plot_gp_model(tracking_pairs[kk][0], tracking_pairs[kk][1],\
					X_DIM,Y_DIM,Z_DIM,figtype_list, use_corrected_cell, self.root_directory)

	def concat_root_to_string(self,path):
		if self.root_directory is not None:
			filename = os.path.join(self.root_directory,path)
			return filename
		else:
			raise Exception('root_directory property must be specified')

	def concat_root_to_all(self,array):
		if self.root_directory is not None:
			filenames = np.array([])
			for item in array:
				filenames = np.append(filenames, os.path.join(self.root_directory,item))
			return filenames
		else:
			raise Exception('root_directory property must be specified')

	def set_inputs(self,filenames_cell,filenames_beads,savefnames,tracking_pairs,fov_dims):
		self.set_filenames_cell(filenames_cell)
		self.set_filenames_beads(filenames_beads)
		self.set_savefnames(savefnames)
		self.set_tracking_pairs(tracking_pairs)
		self.set_fov_dims(fov_dims)

	def set_fov_dims(self,fov_dims):
		self.fov_dims = fov_dims

	def set_filenames_cell(self,filenames_cell):
		self.filenames_cell = filenames_cell

	def get_filenames_cell(self):
		if self.filenames_cell is not None:
			filenames_cell = self.concat_root_to_all(self.filenames_cell)
			dirnames_cell = []
			for item in filenames_cell:
				dirnames_cell.append(os.path.dirname(item))
			return filenames_cell, dirnames_cell 
		else:
			raise Exception('filenames_cell and dirnames_cell must both be specified')
		
	def set_filenames_beads(self,filenames_beads):
		self.filenames_beads = filenames_beads

	def get_filenames_beads(self):
		if self.filenames_beads is not None:
			filenames_beads = self.concat_root_to_all(self.filenames_beads)
			dirnames_beads = []
			for item in filenames_beads:
				dirnames_beads.append(os.path.dirname(item))
			return filenames_beads, dirnames_beads 
		else:
			raise Exception('filenames_beads and dirnames_beads must both be specified')

	def set_out_folder_cell(self,out_folder_cell):
		self.out_folder_cell = out_folder_cell

	def set_out_folder_beads(self,out_folder_beads):
		self.out_folder_beads = out_folder_beads

	def set_savefnames(self,savefnames):
		self.savefnames = savefnames
		
	def get_savenames(self):
		if self.out_folder_cell is None:
			self.out_folder_cell = 'Gel_cell_coords'
		if self.out_folder_beads is None:
			self.out_folder_beads = 'Gel_bead_center_coords'
		out_folder_cell = self.concat_root_to_string(self.out_folder_cell)
		out_folder_beads = self.concat_root_to_string(self.out_folder_beads)
		if not os.path.exists(out_folder_cell):
			os.makedirs(out_folder_cell)
		if not os.path.exists(out_folder_beads):
			os.makedirs(out_folder_beads)
		savefnames_cell = [] 
		savefnames_beads = [] 
		for kk in range(0,len(self.savefnames)):
			savefnames_cell.append(out_folder_cell + '/' + self.savefnames[kk] + '_cell_')
			savefnames_beads.append(out_folder_beads + '/' + self.savefnames[kk] + '_beads.txt') 
		return savefnames_cell, savefnames_beads 

	def set_out_folder(self,out_folder):
		self.out_folder = out_folder

	def set_tracking_pairs(self,tracking_pairs):
		self.tracking_pairs = tracking_pairs

	def get_tracking_pairs(self):
		_ = self.get_out_folder()
		return self.tracking_pairs

	def get_out_folder(self):
		if self.out_folder is None:
			self.out_folder = 'Post_proc_summary'
		out_folder = self.concat_root_to_string(self.out_folder)
		if not os.path.exists(out_folder):
			os.makedirs(out_folder)
		return out_folder

	def get_color_channels(self): #for indexing 3D array with cell image and bead image
		return self.cell_channel, self.bead_channel

	def get_FOV_dims(self):
		fov_dims = self.fov_dims
		X_DIM = fov_dims[0]
		Y_DIM = fov_dims[1]
		Z_DIM = fov_dims[2]
		return X_DIM, Y_DIM, Z_DIM 
		
	def get_cell_thresh(self):
		return self.cell_thresh
		
	def get_tracking_params(self):
		return self.num_feat, self.num_nearest, self.buffer_cell_tracking, self.track_type 

	def get_translation_correction_params(self):
		return self.buffer_cell_translation

	def get_postproc_info(self):
		return self.figtype_list, self.plot_type, self.run_GP, self.use_corrected_cell