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
            String containing the filename format
            Example : input_file='./CytoD/Beads/Gel 2 CytoD%s.tif'
        bead_channel : 
            The color to examine (0=red, 1=green, 2=blue) (our example data uses green for cells)
        X_DIM : float
            Total length of microscope imagery along the x dimension (149.95 for example data)
        Y_DIM : float
            Total length of microscope imagery along the x dimension (149.95 for example data)
        Z_DIM : float
            Total length of microscope imagery along the x dimension (140.0 for example data)

        """

        beads = pre_process.get_bead_centers(filenames_beads,bead_channel, X_DIM, Y_DIM, Z_DIM)
        self.points = beads.points