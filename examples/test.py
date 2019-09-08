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
    (4) running the actual code. As should seem obvious, in this section the actual code
    is run
"""

##########################################################################################
# section (1)
# import the necessary modules to use FM-Track
##########################################################################################

from fmtrack.input_info import input_info
from fmtrack.run_tracking_all_steps import run_tracking_all_steps

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

##########################################################################################
# section (3)
# create an input_info object
##########################################################################################

# in this first line, we initialize the input_info object by passing in the
# root_directory path. The root_directory is the parent folder containing all
# of the subfolders with data inside. This example script runs on the data folder
# included in the examples folder
info = input_info('/Users/<username>/Desktop/data')

# these lines use specific functions to pass in the variables we specified in
# section 2
info.set_filenames_cell(filenames_cell)
info.set_filenames_beads(filenames_beads)
info.set_savefnames(savefnames)
info.set_tracking_pairs(tracking_pairs)

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
