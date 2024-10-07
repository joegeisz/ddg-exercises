// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    int E = mesh->nEdges();
    int V = mesh->nVertices();
    SparseMatrix<size_t> mat(E,V);
    mat.reserve(E*2);
    for(Edge e : mesh->edges()) {
	    int a = e.firstVertex().getIndex();
	    int b = e.secondVertex().getIndex();
	    int e_index = e.getIndex();
	    mat.insert(e_index,a) = 1;
	    mat.insert(e_index,b) = 1;
    }
    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    int F = mesh->nFaces();
    int E = mesh->nEdges();
    SparseMatrix<size_t> mat(F,E);
    mat.reserve(F*3);
    for(Face f : mesh->faces()) {
	    int fi = f.getIndex();
	    for( Edge e : f.adjacentEdges()) {
		    int ei = e.getIndex();
		    mat.insert(fi, ei) = 1;
	    }
    }
    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nVertices());
    for(size_t i: subset.vertices){
	vec(i) = 1;
    }
    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nEdges());
    for(size_t i: subset.edges){
	    vec(i) = 1;
    }
    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nFaces());
    for(size_t i: subset.faces){
	    vec(i) = 1;
    }
    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    MeshSubset star;
    Vector<size_t> V = buildVertexVector(subset);
    Vector<size_t> E = buildEdgeVector(subset);
    Vector<size_t> F = buildFaceVector(subset);
    E = E + this->A0*V;
    F = F + this->A1*E + this->A1*this->A0*V;
    for(size_t i = 0; i < mesh->nVertices(); i++){
	if(V[i] != 0){
	    star.addVertex(i);
	}
    }
    for(size_t i = 0; i < mesh->nEdges(); i++){
	if(E[i] != 0){
	    star.addEdge(i);
	}
    }
    for(size_t i = 0; i < mesh->nFaces(); i++){
	if(F[i] != 0){
	    star.addFace(i);
	}
    }
    return star; 

}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    MeshSubset star;
    Vector<size_t> V = buildVertexVector(subset);
    Vector<size_t> E = buildEdgeVector(subset);
    Vector<size_t> F = buildFaceVector(subset);
    V = V + this->A0.transpose()*E + this->A0.transpose()*this->A1.transpose()*F;
    E = E + this->A1.transpose()*F;
    for(size_t i = 0; i < mesh->nVertices(); i++){
	if(V[i] != 0){
	    star.addVertex(i);
	}
    }
    for(size_t i = 0; i < mesh->nEdges(); i++){
	if(E[i] != 0){
	    star.addEdge(i);
	}
    }
    for(size_t i = 0; i < mesh->nFaces(); i++){
	if(F[i] != 0){
	    star.addFace(i);
	}
    }
    return star; 

}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    MeshSubset st = star(subset);
    MeshSubset cl = closure(subset);
    MeshSubset stcl = star(cl);
    MeshSubset clst = closure(st);
    MeshSubset link(clst.vertices,clst.edges,clst.faces);
    for(size_t v: stcl.vertices){
	link.deleteVertex(v);
    }
    for(size_t e: stcl.edges){
	link.deleteEdge(e);
    }
    for(size_t f: stcl.faces){
	link.deleteFace(f);
    }
    return link;
}


/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    MeshSubset cl = closure(subset);
    if(cl.vertices == subset.vertices && cl.edges == subset.edges && cl.faces == subset.faces){
	return true;
    }
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if(!isComplex(subset)){
	return -1;
    }
    //else if(
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}
