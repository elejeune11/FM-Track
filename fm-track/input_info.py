##########################################################################################
# import packages
##########################################################################################
import os 
import numpy as np
import json

##########################################################################################
# class to get path names and constants that can be modified
##########################################################################################

class input_info():

	def __init__():
		self.assign_defaults()

	def assign_defaults():

		# For get_filenames_cell()
		self.filenames_cell = [ \
			'data/7_16_Magnet/Magnet/Magnetic Bead/Magnet%s.tif',\
			'data/7_16_Magnet/No Magnet/Magnetic Bead/No Magnet%s.tif',\
			'data/7_16_Magnet/No Magnet 2/Magnetic Bead/No Magnet 2%s.tif']
		self.dirnames_cell = [ \
			'data/7_16_Magnet/Magnet/Magnetic Bead/',\
			'data/7_16_Magnet/No Magnet/Magnetic Bead/',\
			'data/7_16_Magnet/No Magnet 2/Magnetic Bead/']

		# For get_filenames_beads()
		self.filenames_beads = [ \
			'data/7_16_Magnet/Magnet/Tracking Beads/Magnet%s.tif',\
			'data/7_16_Magnet/No Magnet/Tracking Beads/No Magnet%s.tif',\
			'data/7_16_Magnet/No Magnet 2/Tracking Beads/No Magnet 2%s.tif']
		self.dirnames_beads = [ \
			'data/7_16_Magnet/Magnet/Tracking Beads/',\
			'data/7_16_Magnet/No Magnet/Tracking Beads/',\
			'data/7_16_Magnet/No Magnet 2/Tracking Beads/']

		# For get_savenames()
		self.out_folder_cell = 'Gel_cell_coords'
		self.out_folder_beads = 'Gel_bead_center_coords'
		self.savefnames = [ \
			'data/7_16_Magnet',\
			'data/7_16_NoMagnet',\
			'data/7_16_NoMagnet2']

		# For get_tracking_pairs()
		self.out_folder = 'Post_proc_summary'
		self.tracking_pairs = [ \
			['data/7_16_Magnet','../7_16_NoMagnet'],\
			['data/7_16_NoMagnet','../7_16_NoMagnet2'],\
			['data/7_16_Magnet','../7_16_NoMagnet2']]

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


	def get_filenames_cell():
		return self.filenames_cell, self.dirnames_cell 
		
	def get_filenames_beads():
		return self.filenames_beads, self.dirnames_beads
		
	def get_savenames():
		if not os.path.exists(self.out_folder_cell):
			os.makedirs(self.out_folder_cell)
		if not os.path.exists(self.out_folder_beads):
			os.makedirs(self.out_folder_beads)
		savefnames_cell = [] 
		savefnames_beads = [] 
		for kk in range(0,len(self.savefnames)):
			savefnames_cell.append(self.out_folder_cell + '/' + self.savefnames[kk] + '_cell_')
			savefnames_beads.append(self.out_folder_beads + '/' + self.savefnames[kk] + '_beads.txt') 
		return savefnames_cell, savefnames_beads 

	def get_tracking_pairs():
		if not os.path.exists(self.out_folder):
			os.makedirs(self.out_folder)
		return self.tracking_pairs

	def get_color_channels(): #for indexing 3D array with cell image and bead image
		return self.cell_channel, self.bead_channel

	def get_FOV_dims(type):
		fov_dims = self.fov_dims[type-1]
		X_DIM = fov_dims[0]
		Y_DIM = fov_dims[1]
		Z_DIM = fov_dims[2]
		return X_DIM, Y_DIM, Z_DIM 
		
	def get_cell_thresh():
		return self.cell_thresh
		
	def get_tracking_params():
		return self.num_feat, self.num_nearest, self.buffer_cell, self.track_type 

	def get_translation_correction_params():
		return self.buffer_cell 

	def get_postproc_info():
		return self.figtype_list, self.plot_type, self.run_GP, self.use_corrected_cell