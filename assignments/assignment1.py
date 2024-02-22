import trimesh
import numpy as np
import open3d as o3d
import time
from trimesh import graph, grouping
from trimesh.geometry import faces_to_edges
import numpy as np
from itertools import zip_longest
from heapq import heappush, heappop

##########################################################################################################
##########################################################################################################

#HELPER FUNCTIONS FOR LOOP SUBDIVISION WITH HALF_EDGE DATA STRUCTURE

def loop_subdivision(vertices, edges, faces, iterations=1):

    # New Half Edges with the new vertices are made
    def getNewHalfEdges(last_face, face, newv1, newv2, newv3, v_to_he, new_edges):
        v1, v2, v3 = face
        loop_over = [[newv1, newv2], [newv2, newv3], [newv3, newv1],
                     [newv1, v2], [v2, newv2], [newv2, newv1], 
                     [newv3, newv2], [newv2, v3], [v3, newv3],
                     [v1,newv1], [newv1, newv3], [newv3, v1]]
        prev = None
        for i in range(len(loop_over)):
            edge1 = o3d.geometry.HalfEdge()
            edge1.vertex_indices = loop_over[i]
            edge1.triangle_index = last_face+(i//3)
            if tuple((loop_over[i][1], loop_over[i][0])) in v_to_he: # twin exists already
                twin_idx = v_to_he[tuple((loop_over[i][1], loop_over[i][0]))] 
                new_edges[twin_idx].twin = len(new_edges)
                edge1.twin = twin_idx
            if i%3!=0:
                prev.next = len(new_edges)
            elif prev:
                prev.next = len(new_edges) - 3
            prev = edge1
            v_to_he[tuple(loop_over[i])] = len(new_edges)
            new_edges.append(edge1)
        if prev:
            prev.next = len(new_edges) - 3
        return new_edges
    
    # Calculates the odd vertices according to the rules of Loop Subdivision
    def calculateOddVertices(v1, v2, v3, v_to_he):
        # print("Calculating odd vertices ...", v1, " ", v2, " ", v3)
        edge_loop = [v1, v2, v3]
        midpoints = []
        for i in range(len(edge_loop)):
            cur_edge = edges[v_to_he[tuple((edge_loop[i],edge_loop[(i+1)%3]))]]
            if cur_edge.twin!=-1:
                opp_face = faces[edges[cur_edge.twin].triangle_index]
                for v in opp_face:
                    if v!=edge_loop[i] and v!=edge_loop[(i+1)%3]:
                        v4 = v
                odd = 0.125*(vertices[edge_loop[(i+2)%3]]+vertices[v4]) + 0.375*(vertices[edge_loop[i]]+vertices[edge_loop[(i+1)%3]])
            else:
                # boundary case
                odd = 0.5*(edge_loop[i]+edge_loop[(i+1)%3])
            midpoints.append(odd.tolist())
        return midpoints

    # MAIN LOOP SUDIVISION
    for i in range(iterations):
        # print("Iteration: ", i)
        v_to_he = {}
        new_edges = []
        new_faces = []
        new_vertices = []

        for i in range(len(edges)):
            edge = edges[i]
            v1,v2 = edge.vertex_indices
            v_to_he[tuple((v1,v2))] = i
            
        # print("Split each edge and update vertex positions")
        for face in faces:
            v1, v2, v3 = face
            newv1, newv2, newv3 = 0,0,0
            # midpoint1 = (0.5*(vertices[v1]+vertices[v2])).tolist() # last_vertex
            # midpoint2 = (0.5*(vertices[v2]+vertices[v3])).tolist() # last_vertex+1
            # midpoint3 = (0.5*(vertices[v3]+vertices[v1])).tolist() # last_vertex+2
            midpoint1, midpoint2, midpoint3 = calculateOddVertices(v1, v2, v3, v_to_he)
            if midpoint1 in new_vertices:
                # midpoint is already processed by the twin edge
                newv1 = new_vertices.index(midpoint1)
            else:
                # midpoint is added to the new vertices
                new_vertices.append(midpoint1)
                newv1 = len(new_vertices)-1

            if midpoint2 in new_vertices:
                newv2 = new_vertices.index(midpoint2)
            else:
                new_vertices.append(midpoint2)
                newv2 = len(new_vertices)-1

            if midpoint3 in new_vertices:
                newv3 = new_vertices.index(midpoint3)
            else:
                new_vertices.append(midpoint3)
                newv3 = len(new_vertices)-1

            last_face = len(new_faces)
            newv1 += len(vertices)
            newv2 += len(vertices)
            newv3 += len(vertices)
            new_edges = getNewHalfEdges(last_face, face, newv1, newv2, newv3, v_to_he, new_edges)
            
            face1 = [newv1, newv2, newv3] # last_face
            face2 = [newv1, v2, newv2] # last_face+1
            face3 = [newv3, newv2, v3] # last_face+2
            face4 = [v1, newv1, newv3] # last_face+3
            new_faces.append(face1)
            new_faces.append(face2)
            new_faces.append(face3)
            new_faces.append(face4)

        # print("Update old vertices")
        for v in range(len(vertices)):
            connected_edges = [edge.vertex_indices for edge in edges if v in edge.vertex_indices]
            n = len(connected_edges)
            if n>3:
                beta = 3.0 / (8.0 * n)
            else:
                beta = 3.0 / 16.0
            new_position = (1 - n * beta) * vertices[v]

            for edge in connected_edges:
                other_vertex = edge[0] if edge[1] == v else edge[1]
                new_position += beta * vertices[other_vertex]

            vertices[v] = new_position

        vertices = o3d.utility.Vector3dVector(np.vstack([vertices, new_vertices]))
        faces = o3d.utility.Vector3iVector(new_faces)
        edges = new_edges
        
    return vertices, edges, faces

def subdivision_loop(mesh, iterations=1):
    """
    Apply Loop subdivision to the input mesh for the specified number of iterations.
    :param mesh: input mesh
    :param iterations: number of iterations
    :return: mesh after subdivision
    """
    halfedge_mesh = o3d.geometry.HalfEdgeTriangleMesh.create_from_triangle_mesh(mesh)
    halfedges = halfedge_mesh.half_edges
    faces = np.asarray(halfedge_mesh.triangles)
    vertices = np.asarray(halfedge_mesh.vertices)
    new_vertices, new_edges, new_faces = loop_subdivision(vertices, halfedges, faces, iterations)

    triangle_mesh = o3d.geometry.TriangleMesh(new_vertices, new_faces)
    return triangle_mesh

##########################################################################################################
##########################################################################################################

# HELPER FUNCTIONS FOR LOOP SUBDIVISION WITH ADJACENCY LISTS

def subdivision_loop_adjacency(mesh, iterations=1):
    
    for i in range(iterations):
      # prepare geometry for the loop subdivision
      vertices, faces = mesh.vertices, mesh.faces # [N_vertices, 3] [N_faces, 3]
      edges, edges_face = faces_to_edges(faces, return_index=True) # [N_edges, 2], [N_edges]
      edges.sort(axis=1)
      unique, inverse = grouping.unique_rows(edges)
      
      # split edges to interior edges and boundary edges
      edge_inter = np.sort(grouping.group_rows(edges, require_count=2), axis=1)
      edge_bound = grouping.group_rows(edges, require_count=1)
      
      # set also the mask for interior edges and boundary edges
      edge_bound_mask = np.zeros(len(edges), dtype=bool)
      edge_bound_mask[edge_bound] = True
      edge_bound_mask = edge_bound_mask[unique]
      edge_inter_mask = ~edge_bound_mask
      
      odd = vertices[edges[unique]].mean(axis=1) # [N_oddvertices, 3]
      
      # connect the odd vertices with even vertices
      # however, the odd vertices need further updates over it's position
      # we therefore complete this step later afterwards.
      
      ###########
      # Step 2: #
      ###########
      # find v0, v1, v2, v3 and each odd vertex
      # v0 and v1 are at the end of the edge where the generated odd vertex on
      # locate the edge first
      e = edges[unique[edge_inter_mask]]
      # locate the endpoints for each edge
      e_v0 = vertices[e][:, 0]
      e_v1 = vertices[e][:, 1]
      
      # v2 and v3 are at the farmost position of the two triangle
      # locate the two triangle face
      edge_pair = np.zeros(len(edges)).astype(int)
      edge_pair[edge_inter[:, 0]] = edge_inter[:, 1]
      edge_pair[edge_inter[:, 1]] = edge_inter[:, 0]
      opposite_face1 = edges_face[unique]
      opposite_face2 = edges_face[edge_pair[unique]]
      # locate the corresponding edge
      e_f0 = faces[opposite_face1[edge_inter_mask]]
      e_f1 = faces[opposite_face2[edge_inter_mask]]
      # locate the vertex index and vertex location
      e_v2_idx = e_f0[~(e_f0[:, :, None] == e[:, None, :]).any(-1)]
      e_v3_idx = e_f1[~(e_f1[:, :, None] == e[:, None, :]).any(-1)]
      e_v2 = vertices[e_v2_idx]
      e_v3 = vertices[e_v3_idx]
      
      # update the odd vertices based the v0, v1, v2, v3, based the following:
      # 3 / 8 * (e_v0 + e_v1) + 1 / 8 * (e_v2 + e_v3)
      odd[edge_inter_mask] = 0.375 * e_v0 + 0.375 * e_v1 + e_v2 / 8.0 + e_v3 / 8.0
      
      ###########
      # Step 3: #
      ###########
      # find vertex neightbors for even vertices and update accordingly
      neighbors = graph.neighbors(edges=edges[unique], max_index=len(vertices))
      # convert list type of array into a fixed-shaped numpy array (set -1 to empties)
      neighbors = np.array(list(zip_longest(*neighbors, fillvalue=-1))).T
      # if the neighbor has -1 index, its point is (0, 0, 0), so that it is not included in the summation of neighbors when calculating the even
      vertices_ = np.vstack([vertices, [0.0, 0.0, 0.0]])
      # number of neighbors
      k = (neighbors + 1).astype(bool).sum(axis=1)
      
      # calculate even vertices for the interior case
      beta = (40.0 - (2.0 * np.cos(2 * np.pi / k) + 3) ** 2) / (64 * k)
      even = (
          beta[:, None] * vertices_[neighbors].sum(1)
          + (1 - k[:, None] * beta[:, None]) * vertices
      )
      
      ############
      # Step 1+: #
      ############
      # complete the subdivision by updating the vertex list and face list
      
      # the new faces with odd vertices
      odd_idx = inverse.reshape((-1, 3)) + len(vertices)
      new_faces = np.column_stack(
          [
              faces[:, 0],
              odd_idx[:, 0],
              odd_idx[:, 2],
              odd_idx[:, 0],
              faces[:, 1],
              odd_idx[:, 1],
              odd_idx[:, 2],
              odd_idx[:, 1],
              faces[:, 2],
              odd_idx[:, 0],
              odd_idx[:, 1],
              odd_idx[:, 2],
          ]
      ).reshape((-1, 3)) # [N_face*4, 3]

      # stack the new even vertices and odd vertices
      new_vertices = np.vstack((even, odd)) # [N_vertex+N_edge, 3]
      mesh = trimesh.Trimesh(new_vertices, new_faces)
    
    return mesh

##########################################################################################################
##########################################################################################################

# HELPER FUNCTIONS FOR DECIMATION 

def homogeneous_coordinates_of_plane(a, b, c, x0, y0, z0):
    d = -(a * x0 + b * y0 + c * z0)
    return  d

def computeK(a,b,c,d):
    K = np.zeros((4,4))
    K[0][0] = a**2
    K[0][1] = K[1][0] = a*b
    K[0][2] = K[2][0] = a*c
    K[0][3] = K[3][0] = a*d
    K[1][1] = b**2
    K[1][2] = K[2][1] = b*c
    K[1][3] = K[3][1] = b*d
    K[2][2] = c**2
    K[2][3] = K[3][2] = c*d
    K[3][3] = d**2
    return K

def initK(mesh, KforVertex):
    # get K for all vertices
    for i in range(len(mesh.faces)):
        a, b, c = mesh.face_normals[i]
        v1, v2, v3 = mesh.vertices[mesh.faces[i][0]]
        d = homogeneous_coordinates_of_plane(a, b, c, v1, v2, v3)
        K = computeK(a,b,c,d)
        v1, v2, v3 = mesh.faces[i]
        KforVertex[v1] += K
        KforVertex[v2] += K
        KforVertex[v3] += K
    return KforVertex
        
def collapseEdge(mesh, v1, v2, x, K, vmask):

    # print("Collapsing v1: ", v1, " v2: ", v2, "x: ",x)
    vertices = np.asarray(mesh.vertices)
    faces = np.asarray(mesh.faces)
    vertices = np.vstack([vertices, x])
    mask = np.ones(len(faces), dtype=bool)
    for i in range(len(faces)):
        if (v1 in faces[i] and v2 in faces[i]):
            mask[i] = False
        elif (v1 in faces[i]):
            idx = np.where(faces[i] == v1)[0][0]
            faces[i][idx] = len(vertices)-1
        elif v2 in faces[i]:
            idx = np.where(faces[i] == v2)[0][0]
            faces[i][idx] = len(vertices)-1
    mask = np.transpose(mask)
    
    mesh.vertices = vertices
    mesh.faces = faces
    mesh.update_faces(mask)
    mesh.remove_unreferenced_vertices()
    
    return mesh

def simplify_quadric_error(mesh, face_count=1):
    """
    Apply quadratic error mesh decimation to the input mesh until the target face count is reached.
    :param mesh: input mesh
    :param face_count: number of faces desired in the resulting mesh.
    :return: mesh after decimation
    """
    while len(mesh.faces)>face_count:
        QForEdge = []
        KforVertex = []
        KForEdge = {}
        for _ in mesh.vertices:
            KforVertex.append(np.zeros((4,4)))
        KforVertex = initK(mesh, KforVertex) #get K for all vertices
        for i in range(len(mesh.edges)):
            v1,v2 = mesh.edges[i]

            K = KforVertex[v1] + KforVertex[v2]
            Q = K.copy()
            Q[3,:] = np.array([0, 0, 0, 1])
            vec = np.zeros((4,1))
            vec[3] = 1
            try:
                Qinv = np.linalg.inv(Q)
                x = np.matmul(Qinv, vec)
            except:
                print("Error: Matrix doesn't have an inverse")
                x = (np.array(mesh.vertices[v1]) + np.array(mesh.vertices[v2]))/2.0
                x = np.vstack([x[:, np.newaxis], [1]])

            cost = np.matmul(np.matmul(np.transpose(x), K), x).reshape((1))
            KForEdge[i] = (K,x)
            heappush(QForEdge,[cost, i]) # x is the new vertex and i is the edge index

        vmask = np.ones(len(mesh.vertices), dtype=bool)
        min_cost, idx = heappop(QForEdge)
        K_foredge, x = KForEdge[idx]
        v1, v2 = mesh.edges[idx]
        mesh = collapseEdge(mesh, v1, v2, x[0:3].reshape((3)), K_foredge, vmask)
        
    if len(mesh.faces)>face_count:
        print("NOT POSSIBLE")
        return mesh
    return mesh


##########################################################################################################
##########################################################################################################

if __name__ == '__main__':
    
    # Loop Subdivision using Half Edge Data Structure

    print("Loop Subdivision using Half Edge Data Structure")
    print('-------------------------------------------------')
    mesh = o3d.io.read_triangle_mesh('assets/cube.obj')
    print(f'Mesh Info: {mesh}')
    print("1 Iteration")
    start = time.time()
    mesh_subdivided = subdivision_loop(mesh, iterations=1)
    end = time.time()
    print(f'Subdivided Mesh Info: {mesh_subdivided}')
    print("Time taken: ", end-start)
    o3d.io.write_triangle_mesh('assets/assignment1/cube_subdivided_half_edge.obj', mesh_subdivided)    
    print('-------------------------------------------------')
    print("2 Iterations")
    start = time.time()
    mesh_subdivided = subdivision_loop(mesh, iterations=2)
    end = time.time()
    print(f'Subdivided Mesh Info: {mesh_subdivided}')
    print("Time taken: ", end-start)
    o3d.io.write_triangle_mesh('assets/assignment1/cube_subdivided_half_edge_2.obj', mesh_subdivided)    
    print('-------------------------------------------------')
    print("3 Iterations")
    start = time.time()
    mesh_subdivided = subdivision_loop(mesh, iterations=3)
    end = time.time()
    print(f'Subdivided Mesh Info: {mesh_subdivided}')
    print("Time taken: ", end-start)
    o3d.io.write_triangle_mesh('assets/assignment1/cube_subdivided_half_edge_3.obj', mesh_subdivided)    
    print('-------------------------------------------------')

    ##########################################################################################################

    print("Loop Subdivision using Adjacency lists")
    print('-------------------------------------------------')
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    print(f'Mesh Info: {mesh}')
    print("1 Iteration")
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    start = time.time()
    mesh_subdivided = subdivision_loop_adjacency(mesh, 1)
    end = time.time()
    print(f'Subdivided Mesh Info: {mesh_subdivided}')
    print("Time taken: ", end-start)
    mesh_subdivided.export('assets/assignment1/cube_subdivided_adjacency.obj')
    print('-------------------------------------------------')
    print("2 Iterations")
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    start = time.time()
    mesh_subdivided = subdivision_loop_adjacency(mesh, 2)
    end = time.time()
    print(f'Subdivided Mesh Info: {mesh_subdivided}')
    print("Time taken: ", end-start)
    mesh_subdivided.export('assets/assignment1/cube_subdivided_adjacency_2.obj')
    print('-------------------------------------------------')
    print("3 Iterations")
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    start = time.time()
    mesh_subdivided = subdivision_loop_adjacency(mesh, 3)
    end = time.time()
    print(f'Subdivided Mesh Info: {mesh_subdivided}')
    print("Time taken: ", end-start)
    mesh_subdivided.export('assets/assignment1/cube_subdivided_adjacency_3.obj')
    print('-------------------------------------------------')

    ##########################################################################################################

    print("Quadric Decimation")
    print('-------------------------------------------------')
    print(" Face Count: 8")
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    start = time.time()
    mesh_decimated = simplify_quadric_error(mesh, face_count=8)
    end = time.time()
    print(f'Decimated Mesh Info: {mesh_decimated}')
    print("Time taken: ", end-start)
    mesh_decimated.export('assets/assignment1/cube_decimated_8faces.obj')
    print('-------------------------------------------------')
    print(" Face Count: 6")
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    start = time.time()
    mesh_decimated = simplify_quadric_error(mesh, face_count=6)
    end = time.time()
    print(f'Decimated Mesh Info: {mesh_decimated}')
    print("Time taken: ", end-start)
    mesh_decimated.export('assets/assignment1/cube_decimated_6faces.obj')
    print('-------------------------------------------------')
    print(" Face Count: 4")
    mesh = trimesh.creation.box(extents=[1, 1, 1])
    start = time.time()
    mesh_decimated = simplify_quadric_error(mesh, face_count=4)
    end = time.time()
    print(f'Decimated Mesh Info: {mesh_decimated}')
    print("Time taken: ", end-start)
    mesh_decimated.export('assets/assignment1/cube_decimated.obj')