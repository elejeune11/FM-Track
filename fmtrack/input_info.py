##########################################################################################
# import packages
##########################################################################################
import os 
import numpy as np

##########################################################################################
# class to get path names and constants that can be modified
##########################################################################################

class input_info:

	def __init__(self,root_directory=None):

		if root_directory is not None:
			self.root_directory = root_directory
		
		self.assign_defaults()

	def assign_defaults(self):

		# For get_color_channels()
		self.cell_channel = 0 #CellBrite Red
		self.bead_channel = 1 #Green fluorescent beads 

		# For get FOV_dims()
		self.fov_dims = np.array([[149.95, 149.95, 140.0], [141.70, 141.70, 120.0], [149.95, 149.95, 120.0]])

		# For get_cell_thresh()
		self.cell_thresh = 1.0

		# For get_tracking_params()
		self.num_feat = 5
		self.num_nearest = 15
		self.buffer_cell = 0
		self.track_type = 2 # type 1 will NOT perform translation correction, type 2 will

		# For get_translation_correction_params()
		self.buffer_cell = 30

		# For get_postproc_info()
		self.figtype_list = ['.png'] 
		self.plot_type = 6.0
		self.run_GP = False
		self.use_corrected_cell = True

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

	def get_FOV_dims(self, type):
		fov_dims = self.fov_dims[type-1]
		X_DIM = fov_dims[0]
		Y_DIM = fov_dims[1]
		Z_DIM = fov_dims[2]
		return X_DIM, Y_DIM, Z_DIM 
		
	def get_cell_thresh(self):
		return self.cell_thresh
		
	def get_tracking_params(self):
		return self.num_feat, self.num_nearest, self.buffer_cell, self.track_type 

	def get_translation_correction_params(self):
		return self.buffer_cell 

	def get_postproc_info(self):
		return self.figtype_list, self.plot_type, self.run_GP, self.use_corrected_cell