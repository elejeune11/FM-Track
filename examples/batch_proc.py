"""
    batch_proc.py: a simplified version of test.py that demonstrates how to
    test a batch of data (multiple data folders)
"""

# import necessary modules
from fmtrack.input_info import input_info
from fmtrack.run_tracking_all_steps import run_tracking_all_steps
import numpy as np
from glob import glob

# set input parameters
filenames_cell = [ \
    'Magnet/Magnetic Bead/Magnet%s.tif',\
    'No Magnet/Magnetic Bead/No Magnet%s.tif',\
    'No Magnet 2/Magnetic Bead/No Magnet 2%s.tif']
filenames_beads = [ \
    'Magnet/Tracking Beads/Magnet%s.tif',\
    'No Magnet/Tracking Beads/No Magnet%s.tif',\
    'No Magnet 2/Tracking Beads/No Magnet 2%s.tif']

savefnames = [ \
    'MagnetSave',\
    'NoMagnetSave',\
    'NoMagnetSave2']

tracking_pairs = [ \
    ['MagnetSave', 'NoMagnetSave'],\
    ['NoMagnetSave', 'NoMagnetSave2'],\
    ['MagnetSave', 'NoMagnetSave2']]

# create the fov_dims variable, which is a 3d array of dimension data. In
# this example, we simply repeat the same numbers three times. Obviously,
# with real data you should modify these numbers to match your experiment
fov_dims1 = np.array([[149.95, 149.95, 140.0], [141.70, 141.70, 120.0], [149.95, 149.95, 120.0]])
fov_dims2 = np.array([[149.95, 149.95, 140.0], [141.70, 141.70, 120.0], [149.95, 149.95, 120.0]])
fov_dims3 = np.array([[149.95, 149.95, 140.0], [141.70, 141.70, 120.0], [149.95, 149.95, 120.0]])
fov_dims = np.dstack((fov_dims1,fov_dims2,fov_dims3))

# list out root_directories from parent folder. The star will be replaced by the name
# of all subdirectories in the <parent-folder>. For example, if you have 10 data folders
# between "data1" and "data10", glob will automatically change the star to equal each
# of these names and return an array of all of them to loop through
root_directories = glob('/Users/<username>/Desktop/<parent-folder>/*')

for i in range(fov_dims.shape[0]):
    print('Running on ' + root_directories[i])
    info = input_info(root_directories[i])
    info.set_inputs(filenames_cell, filenames_beads, savefnames, tracking_pairs, fov_dims[i])
    run_tracking_all_steps(True,True,True,info)
