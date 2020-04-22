import numpy as np 
from pathlib import Path
from . import pre_process

class FMBeads:

    def __init__(self,points=None):
        self.points = points

    # saves bead positions as text file
    def save_as_txt(self,filename):
        filename = Path(filename).with_suffix('.txt')
        np.savetxt(filename,self.points)

    # loads bead positions from text file
    def load_from_txt(self,filename):
        self.points = np.loadtxt(filename)

    # load from x_pos, y_pos, and z_pos
    def load_from_positions(self, x_pos, y_pos, z_pos):
        self.points = np.transpose(np.vstack((x_pos,y_pos,z_pos)))

    # get x_pos, y_pos, and z_pos
    def get_xyz(self):
        return self.points[:,0], self.points[:,1], self.points[:,2]

    def get_bead_centers(self, filenames_beads,bead_channel, X_DIM, Y_DIM, Z_DIM):
        """Creates a FMBeads object from image data

        Parameters
        ----------
        filenames_beads : str
        bead_channel : str
        X_DIM : float
            Total length of microscope imagery along the x dimension
        Y_DIM : float
            Total length of microscope imagery along the x dimension
        Z_DIM : float
            Total length of microscope imagery along the x dimension

        """

        beads = pre_process.get_bead_centers(filenames_beads,bead_channel, X_DIM, Y_DIM, Z_DIM)
        self.points = beads.points