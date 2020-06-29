"""
    simple_example.py: this example file uses the most basic objects from FM-Track
    to demonstrate how it can be used flexibly. Feel free to modify the synthetic deformation
    by changing the displace function below.
"""

# import necessary modules
import fmtrack
import numpy as np 

# this file path is the folder where the data will be saved
tracker_filepath = './example_folder'

# set the filepath to the .msh file using mesh_init_filepath
mesh_init_filepath = ''

##########################################################################################
# section (1)
# build a custom bead dataset
##########################################################################################

# define a displacement function
def displace(pos_orig):
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
    points_final = np.append(points_final, np.array([displace(bead)]), axis=0)

# here we initialize two bead data sets using the FMBeads() class
beads_init = fmtrack.FMBeads(points_init)
beads_final = fmtrack.FMBeads(points_final)


##########################################################################################
# section (2)
# import the mesh files
##########################################################################################

# import the initial mesh
mesh_init = fmtrack.FMMesh()
mesh_init.import_msh_file(mesh_init_filepath)

# use the displace function to simulate where mesh_final would be
points_final = np.zeros(mesh_init.points.shape)
for i in range(mesh_init.points.shape[0]):
    points_final[i] = displace(mesh_init.points[i])

mesh_final = fmtrack.FMMesh(points=points_final, faces=mesh_init.faces)


##########################################################################################
# section (3)
# pass these variables into a tracker object and run the tracking algorithm
##########################################################################################

# initialize a FMTracker() object using the data specified
tracker = fmtrack.FMTracker(cell_init=mesh_init, cell_final=mesh_final, beads_init=beads_init, beads_final=beads_final)

# run the tracking algorithm
tracker.run_tracking()

# save the tracker so the data can be re-used later
tracker.save_all(tracker_filepath)