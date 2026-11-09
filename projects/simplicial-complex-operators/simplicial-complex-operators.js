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
                let result = new MeshSubset();

                for (let i of subset.vertices) {
                        result.addVertex(i);
                        let v = this.mesh.vertices[i];
                        for (let e of v.adjacentEdges()) {
                                result.addEdge(e.index);
                        }
                        for (let f of v.adjacentFaces()) {
                                result.addFace(f.index);
                        }
                }

                for (let i of subset.edges) {
                        result.addEdge(i);
                        let h = this.mesh.edges[i].halfedge;
                        for (let he of [h, h.twin]) {
                                if (!he.onBoundary) {
                                        result.addFace(he.face.index);
                                }
                        }
                }

                for (let i of subset.faces) {
                        result.addFace(i);
                }

                return result;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let result = new MeshSubset();

                for (let i of subset.vertices) {
                        result.addVertex(i);
                }

                for (let i of subset.edges) {
                        result.addEdge(i);
                        let h = this.mesh.edges[i].halfedge;
                        result.addVertex(h.vertex.index);
                        result.addVertex(h.twin.vertex.index);
                }

                for (let i of subset.faces) {
                        result.addFace(i);
                        let f = this.mesh.faces[i];
                        for (let e of f.adjacentEdges()) {
                                result.addEdge(e.index);
                                result.addVertex(e.halfedge.vertex.index);
                                result.addVertex(e.halfedge.twin.vertex.index);
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
                let result = new MeshSubset();

                for (let i of subset.vertices) {
                        let v = this.mesh.vertices[i];
                        for (let u of v.adjacentVertices()) {
                                result.addVertex(u.index);
                        }
                        for (let f of v.adjacentFaces()) {
                                for (let e of f.adjacentEdges()) {
                                        if (e.halfedge.vertex !== v && e.halfedge.twin.vertex !== v) {
                                                result.addEdge(e.index);
                                        }
                                }
                        }
                }

                for (let i of subset.edges) {
                        let h = this.mesh.edges[i].halfedge;
                        for (let he of [h, h.twin]) {
                                if (!he.onBoundary) {
                                        for (let v of he.face.adjacentVertices()) {
                                                if (v !== h.vertex && v !== h.twin.vertex) {
                                                        result.addVertex(v.index);
                                                }
                                        }
                                }
                        }
                }

                let forbidden = new Set(subset.vertices);
                for (let i of subset.edges) {
                        let h = this.mesh.edges[i].halfedge;
                        forbidden.add(h.vertex.index);
                        forbidden.add(h.twin.vertex.index);
                }
                for (let i of subset.faces) {
                        for (let v of this.mesh.faces[i].adjacentVertices()) {
                                forbidden.add(v.index);
                        }
                }

                for (let i of Array.from(result.vertices)) {
                        if (forbidden.has(i)) {
                                result.deleteVertex(i);
                        }
                }
                for (let i of Array.from(result.edges)) {
                        let h = this.mesh.edges[i].halfedge;
                        if (forbidden.has(h.vertex.index) || forbidden.has(h.twin.vertex.index)) {
                                result.deleteEdge(i);
                        }
                }

                return result;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // TODO
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // TODO
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let result = new MeshSubset();

                if (subset.faces.size > 0) {
                        let incidentFaces = new Map();
                        for (let i of subset.faces) {
                                for (let e of this.mesh.faces[i].adjacentEdges()) {
                                        incidentFaces.set(e.index, (incidentFaces.get(e.index) || 0) + 1);
                                }
                        }
                        for (let [i, count] of incidentFaces) {
                                if (count % 2 === 1) {
                                        result.addEdge(i);
                                        let h = this.mesh.edges[i].halfedge;
                                        result.addVertex(h.vertex.index);
                                        result.addVertex(h.twin.vertex.index);
                                }
                        }
                } else if (subset.edges.size > 0) {
                        let degree = new Map();
                        for (let i of subset.edges) {
                                let h = this.mesh.edges[i].halfedge;
                                for (let v of [h.vertex.index, h.twin.vertex.index]) {
                                        degree.set(v, (degree.get(v) || 0) + 1);
                                }
                        }
                        for (let [i, count] of degree) {
                                if (count % 2 === 1) {
                                        result.addVertex(i);
                                }
                        }
                }

                return result;
        }
}
