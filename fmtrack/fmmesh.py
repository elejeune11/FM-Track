import numpy as np
import meshio
import trimesh
from . import pre_process
import os

class FMMesh:

    def __init__(self,points=None,faces=None,normals=None,center=None,vol=None):
        self.points = points
        self.faces = faces
        self.normals = normals

    def import_msh_file(self,filename):
        mesh = meshio.read(filename)
        mesh = trimesh.Trimesh(mesh.points, mesh.get_cells_type('triangle'))
        self.points = np.asarray(mesh.vertices)
        self.faces = np.asarray(mesh.faces)
        self.center = mesh.center_mass
        self.vol = mesh.volume
        self.calculate_normals()

    def save_msh_file(self,filename):
        cells = {'triangle':self.faces}
        mesh = meshio.Mesh(self.points,cells)
        meshio.write(filename,mesh)
	
    def calculate_normals(self):
        if self.faces is not None:
            self.normals = np.zeros(self.faces.shape)
            for i in range(self.faces.shape[0]):
                point1 = self.points[self.faces[i,0],:]
                point2 = self.points[self.faces[i,1],:]
                point3 = self.points[self.faces[i,2],:]

                self.normals[i] = calculate_normal(point1, point2, point3)
                """
                if np.dot(self.normals[i], point1 - self.center) > 0:
                    self.normals[i] = -self.normals[i]
                """
        else:
            raise Exception('FMMesh.faces must be defined before calculating normals')

    def calculate_center(self):
        mesh = trimesh.Trimesh(self.points, self.faces)
        self.center = mesh.center_mass

    def calculate_vol(self):
        mesh = trimesh.Trimesh(self.points, self.faces)
        self.vol = mesh.volume

    def get_cell_surface(self, filenames_cell, dims, cell_channel=0, cell_threshold=1.0):
        """Creates object from image data

        Parameters
        ----------
        filenames_cell : str
            String containing the path and filename format
            Example : './CytoD/Cell/Gel 2 CytoD%s.tif'
        dims : np.array
		    Total length of microscope imagery along the x, y, and z dimensions (149.95, 149.95, and 140.0 for the example data)
        cell_channel : int
            The color to examine (0=red, 1=green, 2=blue) (our example data uses red for cells)
        cell_threshold : float
            Minimum voxel color intensity for consideration as part of the cell (1.0 for example data)

        """
        
        mesh = pre_process.get_cell_surface(filenames_cell, dims, color_idx=cell_channel, cell_threshold=cell_threshold)
        self.points = mesh.points
        self.faces = mesh.faces
        self.normals = mesh.normals
        self.center = mesh.center
        self.vol = mesh.vol

    def import_native_files(self,root):
        self.import_points(root + 'vertices.txt')
        self.import_faces(root + 'faces.txt')
        self.import_normals(root + 'normals.txt')
        self.import_center(root + 'center.txt')
        self.import_vol(root + 'volume.txt')

    def import_points(self,filename):
        self.points = np.loadtxt(filename)

    def import_faces(self,filename):
        self.faces = np.loadtxt(filename)

    def import_normals(self,filename):
        self.normals = np.loadtxt(filename)

    def import_center(self,filename):
        self.center = np.loadtxt(filename)

    def import_vol(self,filename):
        self.vol = np.loadtxt(filename)

    def export_points(self,filename):
        np.savetxt(filename, self.points)

    def export_faces(self,filename):
        np.savetxt(filename, self.faces)

    def export_normals(self,filename):
        if self.normals is None:
            self.calculate_normals()
        np.savetxt(filename, self.normals)

    def export_center(self,filename):
        np.savetxt(filename, self.center)

    def export_vol(self,filename):
        np.savetxt(filename, np.array([self.vol]))

    def save_native_files(self,root):
        os.makedirs(root, exist_ok=True)

        if self.normals is None:
            self.calculate_normals()
        if self.center is None:
            self.calculate_center()
        if self.vol is None:
            self.calculate_vol()

        self.export_points(root + 'vertices.txt')
        self.export_faces(root + 'faces.txt')
        self.export_normals(root + 'normals.txt')
        self.export_center(root + 'center.txt')
        self.export_vol(root + 'volume.txt')

        

def calculate_normal(point1,point2,point3):
    vec1 = point2 - point1
    vec2 = point3 - point1
    cross = np.cross(vec1,vec2)
    return cross / np.linalg.norm(cross)