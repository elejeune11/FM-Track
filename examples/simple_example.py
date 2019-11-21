"""
    simple_example.py: this example file uses the most basic objects from FM-Track
    to demonstrate how it can be used flexibly. These objects are the same objects implemented
    within the InputInfo class. The primary distinction is that InputInfo uses a custom
    FM-Track file system, whereas using FM-Track's more basic objects will allow you
    to customize your own file interactions.

    In this particular example, we recommend you use data/mesh_normal.msh as your mesh_init and
    data/mesh_cytod.msh as your mesh_final. Since we are using synthetic deformation data, the
    bead displacements will not be a function of the cellular deformation, but for the purposes
    of this example it is still nice to see how to use cellular data.
"""

# import necessary modules
import fmtrack
import numpy as np 

# add file paths to import data from
# for this example, it is suggested to use the two .msh files saved in the data folder
mesh_init_filepath = ''
mesh_final_filepath = ''

# this file path is where the FMTracker object will be saved
tracker_filepath = ''

# this folder path is where the native plots will be saved
plots_folderpath = ''


##########################################################################################
# section (1)
# build a custom bead dataset
##########################################################################################

# define a deformation function
def deformation(pos_orig):
    u = pos_orig[0] * 0.1
    v = pos_orig[1] * 0.2
    w = pos_orig[2] * -0.1
    pos_final = pos_orig + np.array([u,v,w])
    return pos_final

# generate random initial bead position data
x = np.random.rand(400) * 150
y = np.random.rand(400) * 150
z = np.random.rand(400) * 150
points_init = np.transpose(np.vstack((x,y,z)))
points_final = np.empty((0,3))
for bead in points_init:
    points_final = np.append(points_final, np.array([deformation(bead)]), axis=0)

# here we initialize two bead data sets using the FMBeads() class
beads_init = fmtrack.FMBeads(points_init)
beads_final = fmtrack.FMBeads(points_final)


##########################################################################################
# section (2)
# import the mesh files
##########################################################################################

# here we use the FMMesh() class to import the mesh data from the filepaths specified at the top
mesh_init = fmtrack.FMMesh()
mesh_init.import_msh_file(mesh_init_filepath)
mesh_final = fmtrack.FMMesh()
mesh_final.import_msh_file(mesh_final_filepath)


##########################################################################################
# section (3)
# pass these variables into a tracker object and run the tracking algorithm
##########################################################################################

# initialize a FMTracker() object using the data specified
tracker = fmtrack.FMTracker(cell_init=mesh_init, cell_final=mesh_final, beads_init=beads_init, beads_final=beads_final)

# run the tracking algorithm
tracker.run_tracking()

# save the tracker so the data can be re-used later
tracker.save(tracker_filepath)

##########################################################################################
# section (4)
# save the output plots from the tracking algorithm and display the 3D ones
##########################################################################################

# create a FMPlot() object by initializing it with the tracker
plotter = fmtrack.FMPlot(tracker)

# save the native plots (2d plots)
plotter.save_native_plots(plots_folderpath)

# show the 3D plot
plotter.plot()