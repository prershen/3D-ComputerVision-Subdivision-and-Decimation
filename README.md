# Loop Subdivision and Decimation

## Loop Subdivision

### Loop Subdivision Using Half Edge Data structure
Loop Subdivision was implemented using Half edge representation. I have used open3D's HalfEdge data structure. In particular, open3d.geometry.HalfEdge and open3d.geometry.HalfEdgeTriangleMesh (only for loading obj file) was used. The halfedge data structure comprises of properties as follows:

next: Index of the next HalfEdge in the same triangle.
triangle_index: Index of the triangle containing this half edge
twin: Index of the twin HalfEdge
vertex_indices: Index of the ordered vertices forming this half edge
Loop Subdivision Algorithm:

Initialize the mesh with its original vertices and triangles.
For each original vertex, compute the updated vertex position based on the Loop subdivision rule.
Create new vertices at the midpoints of original edges.
For each original triangle, subdivide it into four smaller triangles using the updated and new vertices.
Smooth the updated vertices based on their neighborhood to achieve a smoother result. Odd vertices are modified based on the neighboring 4 vertices. Even vertices are modified based on a beta value and its adjacent vertices.
Repeat the process for a specified number of iterations or until the desired level of subdivision is achieved.

##Loop Subdivision using Adjacency lists
Loop Subdivision can also be done using Adjacency lists. This method involves using lists to process the mesh. We will be using three main lists:

Vertex list: Store triples of coordinates ( x,y,z ) tuples of indices
Edge List: Stores tuples of indices of vertices from the vertex list for each of the edges
Face List: Stores tuples of indices of vertices from the vertex list for each of the faces
When compared to the above implementation of Loop Subdivision using Half Edge Data structure, we see that as the number of faces increases and the robustness of the model increases, the runtime of the Half edge implementation increases. This is due to the time taken to prcoess a mesh in the form of a half edge data structure. But depending on the kind of operation performed, the half edge data structure can be very easy to use and faster for operations like finding neighboring triangles. It can be expected that when the robustness of the model increases, especially when expensive operations have to be performed, half edge data structure will outperform adjacency lists.

##Quadric Error Mesh Decimation
Quadric Error Decimation uses iterative contractions of vertex pairs to simplify models and maintains surface error approximations using quadric matrices. This can facilitate much better approximations, both visually and with respect to geometric error.

###Algorithm Implementation:

Compute the Q matrices for all the initial vertices.
Compute the optimal contraction target v¯ for each edge (v1, v2 ). The error = v_T(Q1 +Q2 )v of this target vertex becomes the cost of contracting that pair.
Place all the pairs in a heap keyed on cost with the minimum cost pair at the top.
Iteratively remove the pair (v1, v2 ) of least cost from the heap, contract this pair, and update the costs of all valid pairs involving v1.

## How to use
```shell
git clone this repository
ls ./ # you should see index.html and README.md showup in the terminal
code ./ # open this folder via vscode locally
# open and right click on the index.html
# select "Open With Live Server" to run the code over localhost.
```
