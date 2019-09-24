"""
    test.py: an example demonstrating how to run FM-Track. It runs on the example
    data placed inside of the examples folder of FM-Track.
    
    This file is intended to be:
    (1) a demonstration of how to use FM-Track in its most simplistic form
    (2) a sample script to test whether your installation has worked properly
    (3) a script that can be easily modified to run FM-Track on your own data

    The script has several distinct sections that will be labelled by number. They
    are as follows:
    (1) import statements. These statements tell the Python interpreter where to find
    FM-Track so that the software can be used
    (2) folder location specification. In this section, the data storage locations are
    specified. For most scientists, data may be stored using different conventions or
    file structures depending on the experiment. FM-Track makes it easy to use the file
    structures and naming conventions that you prefer
    (3) creation of an input_info object. If you are not familiar with object oriented
    programming, an "object" makes it easy to bundle many variables and functions together.
    Here, an input_info object is created, and all of the folder names are passed
    into the object for easy storage and recall.
    (4) running the actual code. As should seem obvious, in this section we run the code
    (5) optional parameter changes. This section demonstrates additional software
    modifications that you can make
"""

# before running this program, reset this root_directory variable to be the path 
# of your data folder
root_directory = '/Users/<username>/Desktop/data'

##########################################################################################
# section (1)
# import the necessary modules to use FM-Track
##########################################################################################

from fmtrack.input_info import input_info
from fmtrack.run_tracking_all_steps import run_tracking_all_steps
import numpy as np

##########################################################################################
# section (2)
# initialize all variables that will be passed into input_info
##########################################################################################

# these strings specify the paths (relative to the root_directory, explained below)
# where FM-Track should look for your cellular data and bead data
filenames_cell = [ \
    'Magnet/Magnetic Bead/Magnet%s.tif',\
    'No Magnet/Magnetic Bead/No Magnet%s.tif',\
    'No Magnet 2/Magnetic Bead/No Magnet 2%s.tif']
filenames_beads = [ \
    'Magnet/Tracking Beads/Magnet%s.tif',\
    'No Magnet/Tracking Beads/No Magnet%s.tif',\
    'No Magnet 2/Tracking Beads/No Magnet 2%s.tif']

# these strings specify the filenames you would like to use for saving your data.
# FM-Track will append other words to these names to create full filenames. An example
# here would be MagnetSave_cell_faces.txt. Here, the code uses the name MagnetSave to
# name the specific file
savefnames = [ \
    'MagnetSave',\
    'NoMagnetSave',\
    'NoMagnetSave2']

# these strings specify the folder names of the tracking pairs you would like to analyze
tracking_pairs = [ \
    ['MagnetSave', 'NoMagnetSave'],\
    ['NoMagnetSave', 'NoMagnetSave2'],\
    ['MagnetSave', 'NoMagnetSave2']]

# this line uses an array to store the dimensions of each gel. In this case, we store
# the dimensions of Magnet, No Magnet, then No Magnet 2
fov_dims = np.array([[149.95, 149.95, 140.0], [141.70, 141.70, 120.0], [149.95, 149.95, 120.0]])

##########################################################################################
# section (3)
# create an input_info object
##########################################################################################

# in this first line, we initialize the input_info object by passing in the
# root_directory path. The root_directory is the parent folder containing all
# of the subfolders with data inside. This example script runs on the data folder
# included in the examples folder. Set this value at the top of the script
info = input_info(root_directory)

# these lines use specific functions to pass in the variables we specified in
# section 2
info.set_inputs(filenames_cell, filenames_beads, savefnames, tracking_pairs, fov_dims)

# this line demonstrates how to reset a tunable parameter. See section 5 for a
# lengthier description
info.num_feat = 5

##########################################################################################
# section (4)
# run all steps of the program
##########################################################################################

# by passing in True or False, you tell the program whether or not to run 
# these three steps:
# (1) the pre-processing, which delineats the cell structure and bead locations
# (2) the actual tracking, which calculates the deformation
# (3) the post_processing, which graphs the results
run_tracking_all_steps(True,True,True,info)


##########################################################################################
# section (5)
# additional, optional steps of the script
##########################################################################################

# if you would like to change any of these parameters when running your script, make
# sure to copy paste these changes above where run_tracking_all_steps() is called.
# For example, above I demonstrated how to change the parameter num_feat by setting
# it to the default value

# the following lines are optional. By default, input_info specifies these values,
# for the folders that store the cell's coordinates, the bead coordinates, and the
# post_processing summary data. Functionality is still included to change them if
# you'd like, but these lines are not necessary
out_folder_cell = 'Gel_cell_coords'
out_folder_beads = 'Gel_bead_center_coords'
out_folder = 'Post_proc_summary'
info.set_out_folder_cell(out_folder_cell)
info.set_out_folder_beads(out_folder_beads)
info.set_out_folder(out_folder)

# use the following lines to change the tunable parameters of the script, if you'd like.
# All of these values are the default values, but you may set them to whatever you'd like
info.cell_channel = 0 #CellBrite Red
info.bead_channel = 1 #Green fluorescent beads 
info.cell_thresh = 1.0
info.num_feat = 5
info.num_nearest = 15
info.buffer_cell_tracking = 0
info.track_type = 2 # type 1 will NOT perform translation correction, type 2 will
info.buffer_cell_translation = 30
info.figtype_list = ['.png'] 
info.plot_type = 6.0
info.run_GP = False
info.use_corrected_cell = True
info.should_plot = True # 3D plotting using PyVista