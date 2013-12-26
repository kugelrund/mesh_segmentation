import bpy
import mathutils
import math
import numpy
import scipy.linalg
import scipy.cluster

delta = 0.02
eta = 0.15

class MeshPolygonAdj:
    """Wrapper for a MeshPolygon. Saves a list of adjacent faces"""

    def __init__(self, polygon):
        self.__poly = polygon
        self.adjPolys = []
        
    @property
    def edge_keys(self):
        return self.__poly.edge_keys
    
    @property
    def vertices(self):
        return self.__poly.vertices
    
    @property
    def normal(self):
        return self.__poly.normal


def _distance(face1, face2):
    """Computes the distance between the two adjacent faces face1 and face2"""
    center1 = sum(face1.vertices)/len(face1.vertices)
    center2 = sum(face2.vertices)/len(face2.vertices)
    return (delta * numpy.linalg.norm(center2 - center1) + (1 - delta) * eta * 
           (1 - math.cos(mathutils.Vector.angle(face1.normal, face2.normal))))

def _create_affinity_matrix(mesh):
    """Create the adjacency matrix of the given mesh"""
    
    faces = [MeshPolygonAdj(face) for face in mesh.polygons]
    # create matrix of the distances first
    W = numpy.identity(len(faces))
    
    # progress bar
    bpy.context.window_manager.progress_begin(0,100)
    progress = 0
    step = 1/len(mesh.edge_keys)
    
    # find adjacent faces
    for edge in mesh.edge_keys:
        bpy.context.window_manager.progress_update(progress) 
        j = None # index of possible adjacent face
        for i, face in enumerate(faces):
            if edge in face.edge_keys:
                if not (j is None or faces[j] in face.adjPolys):
                    face.adjPolys.append(faces[j])
                    faces[j].adjPolys.append(face)
                    W[i,j] = _distance(face, faces[j])
                    W[j,i] = W[i,j]
                    break
                else:
                    j = i
        progress += step
        
    # TODO: for each non adjacent pair of faces find shortest 
    #       path of adjacent faces (dijkstra?)
    
    # change distance entries to similarities
    sigma = W.sum()/(len(faces) * len(faces))
    
    den = 2 * sigma ** 2    # this should be used but without distances of non-
    den = 2                 # adjacent faces gives bad results so use 2 for now
    
    for (i,j), value in numpy.ndenumerate(W):
        if not (value == 0 or value == 1):
            W[i,j] = math.exp(-W[i,j]/den)
            
    return W

def segment_mesh(mesh, k, coefficients, action = None):
    """Segments the given mesh into k clusters and performs the given 
    action for each cluster
    
    """
    
    # set coefficients
    global delta
    global eta
    delta, eta = coefficients
    
    # affinity matrix
    W = _create_affinity_matrix(mesh)
    # degree matrix
    Dsqrt = numpy.sqrt(numpy.diag([1/entry for entry in sum(W, 0)]))
    # graph laplacian
    L = numpy.dot(numpy.dot(Dsqrt, W), Dsqrt)
 
    # get eigenvectors
    l,V = scipy.linalg.eigh(L, eigvals = (L.shape[0] - k, L.shape[0] - 1))
    # normalize
    V = scipy.cluster.vq.whiten(V)
    
    # apply kmeans
    cluster_res,_ = scipy.cluster.vq.kmeans(V, k)
    # get identification vector
    idx,_ = scipy.cluster.vq.vq(V, cluster_res)
    
    if action:
        action(mesh, k, idx)
        
    bpy.context.window_manager.progress_end()
