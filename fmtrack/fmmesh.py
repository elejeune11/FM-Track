import numpy as np
import meshio
import trimesh
from . import pre_process

class FMMesh:

    def __init__(self,points=None,faces=None,normals=None,center=None,vol=None):
        self.points = points
        self.faces = faces
        self.normals = normals
        self.center = center
        self.vol = vol

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

    def get_cell_surface(self, dirnames_cell,filenames_cell,cell_channel, X_DIM, Y_DIM, Z_DIM, cell_threshold):
        mesh = pre_process.get_cell_surface(dirnames_cell,filenames_cell,cell_channel, X_DIM, Y_DIM, Z_DIM, cell_threshold)
        self.points = mesh.points
        self.faces = mesh.faces
        self.normals = mesh.normals
        self.center = mesh.center
        self.vol = mesh.vol

    def import_native_files(self,root):
        self.import_points(root + '_cell_mesh.txt')
        self.import_faces(root + '_cell_faces.txt')
        self.import_normals(root + '_cell_normals.txt')
        self.import_center(root + '_cell_center.txt')
        self.import_vol(root + '_cell_volume.txt')

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
        self.export_points(root + '_cell_mesh.txt')
        self.export_faces(root + '_cell_faces.txt')
        self.export_normals(root + '_cell_normals.txt')
        self.export_center(root + '_cell_center.txt')
        self.export_vol(root + '_cell_volume.txt')

def calculate_normal(point1,point2,point3):
    vec1 = point2 - point1
    vec2 = point3 - point1
    cross = np.cross(vec1,vec2)
    return cross / np.linalg.norm(cross)