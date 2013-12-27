import bpy
import mathutils
import math
import numpy
import scipy.linalg
import scipy.cluster
import scipy.sparse
import scipy.sparse.csgraph

delta = 0.02
eta = 0.15

def _geodesic_distance(face1, face2):
    """Computes the distance between the two adjacent faces face1 and face2"""
    center1 = sum(face1.vertices)/len(face1.vertices)
    center2 = sum(face2.vertices)/len(face2.vertices)
    return numpy.linalg.norm(center2 - center1)


def _angular_distance(face1, face2):
    return eta * (1 - math.cos(mathutils.Vector.angle(face1.normal, 
                                                      face2.normal)))


def _create_affinity_matrix(mesh):
    """Create the adjacency matrix of the given mesh"""
    
    faces = mesh.polygons
    l = len(faces)
    
    # matrix of geodesic distances
    G = scipy.sparse.csr_matrix((l, l), dtype=float)
    maxG = 0
    # matrix of angular distances
    A = scipy.sparse.csr_matrix((l, l), dtype=float)
    maxA = 0
    
    # progress bar
    bpy.context.window_manager.progress_begin(0, 100)
    progress = 0
    step = 1/len(mesh.edge_keys)
    
    # find adjacent faces
    for edge in mesh.edge_keys:
        bpy.context.window_manager.progress_update(progress) 
        j = None # index of possible adjacent face
        for i, face in enumerate(faces):
            if edge in face.edge_keys:
                if not (j is None or G[i,j] != 0):
                    G[i,j] = _geodesic_distance(face, faces[j])
                    A[i,j] = _angular_distance(face, faces[j])
                    G[j,i] = G[i,j]
                    A[i,j] = A[i,j]
                    if G[i,j] > maxG:
                        maxG = G[i,j]
                    if A[i,j] > maxA:
                        maxA = A[i,j]
                    break
                else:
                    j = i
        progress += step
    
    # weight with delta and maximum value
    G.dot(delta/maxG)
    A.dot((1 - delta)/maxA)
    
    # for each non adjacent pair of faces find shortest path of adjacent faces 
    W = scipy.sparse.csgraph.dijkstra(G + A, directed=False)
    
    # change distance entries to similarities
    sigma = W.sum()/(l * l)
    den = 2 * sigma ** 2
    
    for (i,j), value in numpy.ndenumerate(W):
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
