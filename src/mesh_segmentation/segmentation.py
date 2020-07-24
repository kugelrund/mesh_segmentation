import bpy
import mathutils
import math
import numpy
import scipy.linalg
import scipy.cluster
import scipy.sparse
import scipy.sparse.csgraph
import scipy.sparse.linalg

# Controls weight of geodesic to angular distance. Values closer to 0 give
# the angular distance more importance, values closer to 1 give the geodesic
# distance more importance
delta = None

# Weight of convexity. Values close to zero give more importance to concave
# angles, values close to 1 treat convex and concave angles more equally
eta = None


def _face_center(mesh, face):
    """Computes the coordinates of the center of the given face"""
    center = mathutils.Vector()
    for vert in face.vertices:
        center += mesh.vertices[vert].co
    return center/len(face.vertices)


def _geodesic_distance(mesh, face1, face2, edge):
    """Computes the geodesic distance over the given edge between
    the two adjacent faces face1 and face2"""
    edge_center = (mesh.vertices[edge[0]].co + mesh.vertices[edge[1]].co)/2
    return (edge_center - _face_center(mesh, face1)).length + \
           (edge_center - _face_center(mesh, face2)).length


def _angular_distance(mesh, face1, face2):
    """Computes the angular distance of the given adjacent faces"""
    angular_distance = (1 - math.cos(mathutils.Vector.angle(face1.normal,
                                                            face2.normal)))
    if (face1.normal.dot(_face_center(mesh, face2) -
                         _face_center(mesh, face1))) < 0:
        # convex angles are not that bad so scale down distance a bit
        angular_distance *= eta
    return angular_distance


def _create_distance_matrices(mesh):
    """Creates the matrices of the angular and geodesic distances
    between all adjacent faces. The i,j-th entry of the returned
    matrices contains the distance between the i-th and j-th face.
    """

    faces = mesh.polygons
    l = len(faces)

    # number of pairs of adjacent faces
    num_adj = 0

    # map from edge-key to adjacent faces
    adj_faces_map = {}
    # find adjacent faces by iterating edges
    for index, face in enumerate(faces):
        for edge in face.edge_keys:
            if edge in adj_faces_map:
                adj_faces_map[edge].append(index)
            else:
                adj_faces_map[edge] = [index]

    # average G and cumulated A
    avgG = 0
    avgA = 0
    # helping vectors to create sparse matrix later on
    Arow = []
    Acol = []
    Aval = []
    Grow = []
    Gcol = []
    Gval = []
    # iterate adjacent faces and calculate distances
    for edge, adj_faces in adj_faces_map.items():
        if len(adj_faces) == 2:
            i = adj_faces[0]
            j = adj_faces[1]

            Gtemp = _geodesic_distance(mesh, faces[i], faces[j], edge)
            Atemp = _angular_distance(mesh, faces[i], faces[j])
            Gval.append(Gtemp)
            Grow.append(i)
            Gcol.append(j)
            Gval.append(Gtemp)  # add symmetric entry
            Grow.append(j)
            Gcol.append(i)
            Aval.append(Atemp)
            Arow.append(i)
            Acol.append(j)
            Aval.append(Atemp)  # add symmetric entry
            Arow.append(j)
            Acol.append(i)

            avgG += Gtemp
            avgA += Atemp
            num_adj += 1

        elif len(adj_faces) > 2:
            print("Edge with more than 2 adjacent faces: " + str(adj_faces) + "!")

    # create sparse matrices
    # matrix of geodesic distances
    G = scipy.sparse.csr_matrix((Gval, (Grow, Gcol)), shape=(l, l))
    # matrix of angular distances
    A = scipy.sparse.csr_matrix((Aval, (Arow, Acol)), shape=(l, l))

    avgG /= num_adj
    avgA /= num_adj

    # weight with delta and average value
    G = G.dot(delta/avgG)
    A = A.dot((1 - delta)/avgA)

    return G, A


def _create_affinity_matrix(mesh):
    """Create the adjacency matrix of the given mesh"""

    l = len(mesh.polygons)
    print("mesh_segmentation: Creating distance matrices...")
    G, A = _create_distance_matrices(mesh)

    print("mesh_segmentation: Finding shortest paths between all faces...")
    # for each non adjacent pair of faces find shortest path of adjacent faces
    W = scipy.sparse.csgraph.dijkstra(G + A, directed = False)
    inf_indices = numpy.where(numpy.isinf(W))
    W[inf_indices] = 0

    print("mesh_segmentation: Creating affinity matrix...")
    # change distance entries to similarities
    sigma = W.sum()/(l ** 2)
    den = 2 * (sigma ** 2)
    W = numpy.exp(-W/den)
    W[inf_indices] = 0
    numpy.fill_diagonal(W, 1)

    return W


def _initial_guess(Q, k):
    """Computes an initial guess for the cluster-centers"""
    n = Q.shape[0]
    min_value = 2
    min_indices=(-1,-1)
    for (i,j), value in numpy.ndenumerate(Q):
        if i != j and value < min_value:
            min_value = Q[i,j]
            min_indices = (i,j)

    chosen = [min_indices[0], min_indices[1]]
    for _ in range(2,k):
        min_max = float("inf")
        cur_max = 0
        new_index = -1
        for i in range(n):
            if i not in chosen:
                cur_max = Q[chosen,i].max()
                if cur_max < min_max:
                    min_max = cur_max
                    new_index = i
        chosen.append(new_index)

    return chosen


def segment_mesh(mesh, k, coefficients, action, ev_method):
    """Segments the given mesh into k clusters and performs the given
    action for each cluster
    """

    # set coefficients
    global delta
    global eta
    delta, eta = coefficients

    # affinity matrix
    W = _create_affinity_matrix(mesh)
    print("mesh_segmentation: Calculating graph laplacian...")
    # degree matrix
    Dsqrt = numpy.sqrt(numpy.reciprocal(W.sum(1)))
    # graph laplacian
    L = ((W * Dsqrt).transpose() * Dsqrt).transpose()

    print("mesh_segmentation: Calculating eigenvectors...")
    # get eigenvectors
    if ev_method == 'dense':
        _, V = scipy.linalg.eigh(L, eigvals = (L.shape[0] - k, L.shape[0] - 1))
    else:
        _, V = scipy.sparse.linalg.eigsh(L, k)
    # normalize each column to unit length
    V = V / [numpy.linalg.norm(column) for column in V.transpose()]

    print("mesh_segmentation: Preparing kmeans...")
    # compute association matrix
    Q = V.dot(V.transpose())
    # compute initial guess for clustering
    initial_clusters = _initial_guess(Q, k)

    print("mesh_segmentation: Applying kmeans...")
    # apply kmeans
    cluster_res,_ = scipy.cluster.vq.kmeans(V, V[initial_clusters,:])
    # get identification vector
    idx,_ = scipy.cluster.vq.vq(V, cluster_res)

    print("mesh_segmentation: Done clustering!")
    # perform action with the clustering result
    if action:
        action(mesh, k, idx)
