root_directory = '/Users/<username>/Desktop/data'

# section (1)
# import the necessary modules to use FM-Track

from fmtrack.inputinfo import InputInfo
import numpy as np

# section (2)
# initialize all variables that will be passed into InputInfo

filenames_cell = [ \
    'CytoD/Cell/Gel 2 CytoD%s.tif',\
    'Endo1/Cell/Gel 2 Endo1%s.tif',\
    'Normal/Cell/Gel 2 Normal%s.tif']
filenames_beads = [ \
    'CytoD/Beads/Gel 2 CytoD%s.tif',\
    'Endo1/Beads/Gel 2 Endo1%s.tif',\
    'Normal/Beads/Gel 2 Normal%s.tif']

savefnames = [ \
    'CytoDSave',\
    'Endo1Save',\
    'NormalSave']

tracking_pairs = [ \
    ['CytoDSave', 'Endo1Save'],\
    ['CytoDSave', 'NormalSave'],\
    ['Endo1Save', 'NormalSave']]

fov_dims = np.array([149.95, 149.95, 140.0])

# section (3)
# create an InputInfo object

info = InputInfo(root_directory)
info.set_inputs(filenames_cell, filenames_beads, savefnames, tracking_pairs, fov_dims)

# section (4)
# additional, optional steps of the script

#out_folder_cell = 'Gel_cell_coords'
#out_folder_beads = 'Gel_bead_center_coords'
#out_folder = 'Post_proc_summary'
#info.set_out_folder_cell(out_folder_cell)
#info.set_out_folder_beads(out_folder_beads)
#info.set_out_folder(out_folder)

#info.cell_channel = 0 #CellBrite Red
#info.bead_channel = 1 #Green fluorescent beads 
#info.cell_thresh = 1.0
#info.num_feat = 5
#info.num_nearest = 15
#info.buffer_cell_tracking = 0
#info.track_type = 2 # type 1 will NOT perform translation correction, type 2 will
#info.buffer_cell_translation = 30
#info.figtype_list = ['.png'] 
#info.plot_type = 6.0
#info.run_GP = False
#info.use_corrected_cell = True
#info.should_plot = True # toggles 3D plotting using PyVista
#info.print_progress = True # toggles all printing in algorithm

# section (5)
# run all steps of the program

info.run_tracking_all_steps(True,True,True)