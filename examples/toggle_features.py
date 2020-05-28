"""
    toggle_features.py: an example demonstrating how to run FM-Track. It runs on the example
    data placed inside of the examples folder of FM-Track. This file demonstrates
    how to modify many important parameters in FM-Track.
"""


import fmtrack
import numpy as np

should_plot = False   # toggles on and off 3D plotting


# set microscope field-of-view dimensions
dims = np.array([149.95,149.95,140])

# (1) compute cellular boundary from images
print('Importing cell data')
cell_channel = 0
cell_threshold = 1
cell_init = fmtrack.FMMesh()
cell_init.get_cell_surface('./data/CytoD/Cell/Gel 2 CytoD%s.tif', dims, cell_channel=cell_channel, cell_threshold=cell_threshold)
cell_final = fmtrack.FMMesh()
cell_final.get_cell_surface('./data/Normal/Cell/Gel 2 Normal%s.tif', dims, cell_channel=cell_channel, cell_threshold=cell_threshold)

# (2) find bead positions from images
print('Importing bead data')
bead_channel = 1
beads_init = fmtrack.FMBeads()
beads_init.get_bead_centers('./data/CytoD/Beads/Gel 2 CytoD%s.tif', dims, bead_channel=bead_channel)
beads_final = fmtrack.FMBeads()
beads_final.get_bead_centers('./data/Normal/Beads/Gel 2 Normal%s.tif', dims, bead_channel=bead_channel)

# (3) run tracking algorithm
tracker = fmtrack.FMTracker(cell_init=cell_init, cell_final=cell_final, beads_init=beads_init, beads_final=beads_final)
tracker.print_progress = True   # toggles on and off printing
tracker.track_type = 2   # 1 = no tracking correction, 2 = tracking correctiion
tracker.num_feat = 5   # number of feature vectors
tracker.num_nearest = 15   # number of beads in final state to compare to a particular bead in the initial state
tracker.buffer_cell = 60   # buffer around cell beyond which beads are considered "far field"
tracker.use_box = False   # True = far field boundary is rectangular prism, False = far field boundary is determined by dist. to cell surface
tracker.should_remove_spurious = True   # if true, removes spurious far-field beads after tracking. Only use if track_type=2
tracker.run_tracking()

# (4) save all of the output data
tracker.run_gp = False   # if true, makes and plots a gp model 
tracker.save_native_mesh = True   # true = saves cell boundaries as text files, false = saves gmsh files
tracker.save_all('./data/Track_CytoD_to_Normal')