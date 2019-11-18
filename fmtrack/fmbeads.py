import numpy as np 
from pathlib import Path
from . import pre_process

class FMBeads:

    def __init__(self,points=None):
        self.points = points

    def save_as_txt(self,filename):
        filename = Path(filename).with_suffix('.txt')
        np.savetxt(filename,self.points)

    def load_from_txt(self,filename):
        self.points = np.loadtxt(filename)

    def get_bead_centers(self, dirnames_beads,filenames_beads,bead_channel, X_DIM, Y_DIM, Z_DIM):
        self.points = pre_process.get_bead_centers(dirnames_beads,filenames_beads,bead_channel, X_DIM, Y_DIM, Z_DIM)