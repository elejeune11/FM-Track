##########################################################################################
# import packages
##########################################################################################
import os 
import numpy as np
from . import pre_process
from . import post_process
from . import tracking
from . import fmtracker
from . import fmplot
from . import fmbeads
from . import fmmesh

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
		self.run_gp = False
		self.use_corrected_cell = True
		self.should_plot = False
		self.print_progress = True

		self._has_been_constructed = False

		self.set_out_folder_cell('Gel_cell_coords')
		self.set_out_folder_beads('Gel_bead_center_coords')
		self.set_out_folder('Post_proc_summary')

	def run_tracking_all_steps(self,run_pre_process, run_tracking, run_post_process):
		if not run_tracking:
			self.load_pre_processed_data()
			if not run_pre_process:
				self.load_tracking_data()

		if not self._has_been_constructed:
			self._construct_arrays()
		if run_pre_process:
			self.run_pre_process()
		if run_tracking:
			if not run_pre_process:
				self.load_pre_processed_data()
			self.run_tracking()
		if run_post_process:
			self.run_post_process()

	def _construct_arrays(self):
		"""Constructs the following:

			self.pairs_array : dict
				Indexed by the folder name of the tracking pair in the format 'Track_<beads1>_to_<beads2>'. Each element
				contains an array of the form ['<beads1>', '<beads2>'], one row of self.tracking_pairs
			self.tracker_dict : dict
				Indexed in the same way as self.pairs_array. Each element is an FMTracker object
			self.beads_dict : dict
				Indexed by values of self.savefnames. Each element is an FMBeads object containing the bead data for that state
			self.mesh_dict : dict
				Indexed by values of self.savefnames. Each element is an FMMesh object containing the cell data for that state

		"""

		self._has_been_constructed = True
		self.pairs_array = {}
		self.tracker_dict = {}
		for kk in range(0,len(self.tracking_pairs)):
			idx_init = np.where(self.savefnames==self.tracking_pairs[kk][0])[0][0]
			idx_final = np.where(self.savefnames==self.tracking_pairs[kk][1])[0][0]
			filename_pair = foldername = os.path.join(self.root_directory,'Track_'+self.savefnames[idx_init]+'_to_'+self.savefnames[idx_final])
			self.pairs_array[filename_pair] = self.tracking_pairs[kk]
			self.tracker_dict[filename_pair] = fmtracker.FMTracker()

		self.beads_dict = {}
		self.mesh_dict = {}
		for kk in range(0,len(self.savefnames)):
			self.beads_dict[self.savefnames[kk]] = fmbeads.FMBeads()
			self.mesh_dict[self.savefnames[kk]] = fmmesh.FMMesh()


	def run_pre_process(self):
		if not self._has_been_constructed:
			self._construct_arrays()
		if self.print_progress:
			print('Running Pre-Processing Step')

		filenames_beads, dirnames_beads = self.get_filenames_beads()
		filenames_cell, dirnames_cell  = self.get_filenames_cell()
		savefnames_cell, savefnames_beads  = self.get_savenames()
		X_DIM, Y_DIM, Z_DIM = self.get_FOV_dims()

		for kk in range(0,len(filenames_beads)):
			self.beads_dict[self.savefnames[kk]].get_bead_centers(filenames_beads[kk], self.fov_dims, bead_channel=self.bead_channel)
			self.beads_dict[self.savefnames[kk]].save_as_txt(savefnames_beads[kk])
			self.mesh_dict[self.savefnames[kk]].get_cell_surface(filenames_cell[kk], self.fov_dims, cell_channel=self.cell_channel, cell_threshold=self.cell_thresh)
			self.mesh_dict[self.savefnames[kk]].save_native_files(savefnames_cell[kk])

	def load_pre_processed_data(self):
		if not self._has_been_constructed:
			self._construct_arrays()

		savefnames_cell, savefnames_beads  = self.get_savenames()

		for kk in range(0,len(self.savefnames)):
			self.beads_dict[self.savefnames[kk]].load_from_txt(savefnames_beads[kk])
			self.mesh_dict[self.savefnames[kk]].import_native_files(savefnames_cell[kk])
		

	def run_tracking(self):
		if not self._has_been_constructed:
			self._construct_arrays()

		for key in self.pairs_array:
			idx_init = self.pairs_array[key][0]
			idx_final = self.pairs_array[key][1]

			mesh_init = self.mesh_dict[idx_init]
			mesh_final = self.mesh_dict[idx_final]
			bead_init = self.beads_dict[idx_init]
			bead_final = self.beads_dict[idx_final]

			tracker = fmtracker.FMTracker(mesh_init, mesh_final, bead_init, bead_final)
			tracker.num_feat = self.num_feat
			tracker.num_nearest = self.num_nearest
			tracker.buffer_cell = self.buffer_cell_tracking
			tracker.track_type = self.track_type
			tracker.run_tracking()

			foldername = os.path.join(self.root_directory,key)
			tracker.save_res(foldername,True)
			
			if self.run_gp:
				tracker.create_gp_model()
				tracker.save_gp_model(foldername)

			self.tracker_dict[key] = tracker

	def load_tracking_data(self):
		if not self._has_been_constructed:
			self._construct_arrays()

		for key in self.pairs_array:
			tracker = fmtracker.FMTracker()
			foldername = os.path.join(self.root_directory,key)
			cell_init_foldername = os.path.join(self.root_directory,'Gel_cell_coords',self.pairs_array[key][0]+'_cell_')
			cell_final_foldername = os.path.join(self.root_directory,'Gel_cell_coords',self.pairs_array[key][1]+'_cell_')
			tracker.load_res(foldername,cell_init_foldername,cell_final_foldername,True)

			if self.run_gp:
				tracker.load_gp_model(foldername)

			self.tracker_dict[key] = tracker

	def run_post_process(self):
		if not self._has_been_constructed:
			self._construct_arrays()
		if self.print_progress:
			print('Running Post-Processing Step')

		for key in self.pairs_array:
			plotter = fmplot.FMPlot(self.tracker_dict[key])
			plotter.figtype_list = self.figtype_list
			foldername = os.path.join(self.root_directory,key)
			plotter.save_native_plots(foldername)
			if self.should_plot:
				plotter.plot()
			filename = os.path.join(foldername,'translation.png')
			self.tracker_dict[key].save_mars_figure(filename)

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
		self.savefnames = np.asarray(savefnames)
		
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

	def get_FOV_dims(self):
		fov_dims = self.fov_dims
		X_DIM = fov_dims[0]
		Y_DIM = fov_dims[1]
		Z_DIM = fov_dims[2]
		return X_DIM, Y_DIM, Z_DIM 