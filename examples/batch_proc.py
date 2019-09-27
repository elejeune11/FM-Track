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
print(root_directories)

for i in range(fov_dims.shape[0]):
    info = input_info(root_directories[i])
    info.set_inputs(filenames_cell, filenames_beads, savefnames, tracking_pairs, fov_dims[i])
    run_tracking_all_steps(True,True,True,info)