"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                // TOeDO 
                let i = 0;
                for(const v of mesh.vertices) {
                        v.index = i++;
                }
                i = 0;
                for(const e of mesh.edges) {
                        e.index = i++;
                }
                i = 0;
                for(const f of mesh.faces) {
                        f.index = i++;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                //TOeDO

                let T = new Triplet(mesh.edges.length, mesh.vertices.length);
                for (let e of mesh.edges) {
                        T.addEntry(1, e.index, e.halfedge.vertex.index);
                        T.addEntry(1, e.index, e.halfedge.twin.vertex.index);
                }
                return SparseMatrix.fromTriplet(T);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                // TOeDO
                let T = new Triplet(mesh.faces.length, mesh.edges.length);
                for (let f of mesh.faces) {
                        let h = f.halfedge;
                        do {
                                T.addEntry(1, f.index, h.edge.index);
                                h = h.next;
                        } while (h !== f.halfedge);
                }
                return SparseMatrix.fromTriplet(T);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let v = DenseMatrix.zeros(this.mesh.vertices.length, 1);
                for (let i of subset.vertices) {
                        v.set(1, i, 0);
                }
                return v;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let v = DenseMatrix.zeros(this.mesh.edges.length, 1);
                for (let i of subset.edges) {
                        v.set(1, i, 0);
                }
                return v;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let v = DenseMatrix.zeros(this.mesh.faces.length, 1);
                for (let i of subset.faces) {
                        v.set(1, i, 0);
                }
                return v;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                let vertices = this.buildVertexVector(subset);
                let edges = this.buildEdgeVector(subset);
                let faces = this.buildFaceVector(subset);

                let incidentEdges = this.A0.timesDense(vertices);
                let incidentFaces = this.A1.timesDense(edges).plus(this.A1.timesDense(incidentEdges));

                let starEdges = edges.plus(incidentEdges);
                let starFaces = faces.plus(incidentFaces);

                let result = new MeshSubset();
                for (let i = 0; i < vertices.nRows(); i++) {
                        if (vertices.get(i, 0) > 0) {
                                result.addVertex(i);
                        }
                }
                for (let i = 0; i < starEdges.nRows(); i++) {
                        if (starEdges.get(i, 0) > 0) {
                                result.addEdge(i);
                        }
                }
                for (let i = 0; i < starFaces.nRows(); i++) {
                        if (starFaces.get(i, 0) > 0) {
                                result.addFace(i);
                        }
                }

                return result;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let vertices = this.buildVertexVector(subset);
                let edges = this.buildEdgeVector(subset);
                let faces = this.buildFaceVector(subset);

                let edgesOfFaces = this.A1.transpose().timesDense(faces);
                let verticesOfEdges = this.A0.transpose().timesDense(edges);
                let verticesOfFaces = this.A0.transpose().timesDense(edgesOfFaces);

                let closureEdges = edges.plus(edgesOfFaces);
                let closureVertices = vertices.plus(verticesOfEdges).plus(verticesOfFaces);

                let result = new MeshSubset();
                for (let i = 0; i < closureVertices.nRows(); i++) {
                        if (closureVertices.get(i, 0) > 0) {
                                result.addVertex(i);
                        }
                }
                for (let i = 0; i < closureEdges.nRows(); i++) {
                        if (closureEdges.get(i, 0) > 0) {
                                result.addEdge(i);
                        }
                }
                for (let i = 0; i < faces.nRows(); i++) {
                        if (faces.get(i, 0) > 0) {
                                result.addFace(i);
                        }
                }

                return result;
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let closedStar = this.closure(this.star(subset));
                let starClosure = this.star(this.closure(subset));

                for (let i of starClosure.vertices) {
                        closedStar.deleteVertex(i);
                }
                for (let i of starClosure.edges) {
                        closedStar.deleteEdge(i);
                }
                for (let i of starClosure.faces) {
                        closedStar.deleteFace(i);
                }

                return closedStar;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                let closure = this.closure(subset);
                return closure.vertices.size === subset.vertices.size
                        && closure.edges.size === subset.edges.size
                        && closure.faces.size === subset.faces.size;
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                if (!this.isComplex(subset)) {
                        return -1;
                }
                if (subset.faces.size === 0 && subset.edges.size === 0) {
                        return subset.vertices.size > 0 ? 0 : -1;
                }

                let vertices = this.buildVertexVector(subset);
                let edges = this.buildEdgeVector(subset);
                let faces = this.buildFaceVector(subset);

                let degree = subset.faces.size > 0 ? 2 : 1;
                let coveredEdges = degree === 2 ? this.A1.transpose().timesDense(faces) : edges;
                let coveredVertices = this.A0.transpose().timesDense(coveredEdges);

                for (let i = 0; i < edges.nRows(); i++) {
                        if (edges.get(i, 0) > 0 && coveredEdges.get(i, 0) === 0) {
                                return -1;
                        }
                }
                for (let i = 0; i < vertices.nRows(); i++) {
                        if (vertices.get(i, 0) > 0 && coveredVertices.get(i, 0) === 0) {
                                return -1;
                        }
                }

                return degree;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let degree = this.isPureComplex(subset);
                let result = new MeshSubset();

                if (degree === 2) {
                        let faceCounts = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                        for (let i = 0; i < faceCounts.nRows(); i++) {
                                if (faceCounts.get(i, 0) % 2 === 1) {
                                        result.addEdge(i);
                                }
                        }
                        result = this.closure(result);
                } else if (degree === 1) {
                        let degrees = this.A0.transpose().timesDense(this.buildEdgeVector(subset));
                        for (let i = 0; i < degrees.nRows(); i++) {
                                if (degrees.get(i, 0) % 2 === 1) {
                                        result.addVertex(i);
                                }
                        }
                }

                return result;
        }
}
