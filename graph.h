// graph.h <Starter Code>
// Arthor: Hua Xu, UIC, spring 2021
//
// Basic graph class using adjacency list representation. no
// limitions to any nodes;
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//

#pragma once
#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include<map>
#include <unordered_map>
using namespace std;

template<typename VertexT, typename WeightT>
class graph {
private:
    // to store the adjacent node in an unordered_map
    unordered_map<VertexT, map<VertexT, WeightT>> adjList;
    // store the vertex in a set for quickly search
    set<VertexT>  Vertices;
    // store in a vector
    vector<VertexT> vertexs;

public:
    //
    // constructor:
    // NOTE: the graph is implemented using an adjacency list.
    //
    graph() {
        adjList.clear();
        Vertices.clear();
    }
    //
    // NumVertices
    //
    // Returns the # of vertices currently in the graph.
    //
    int NumVertices() const {
        return static_cast<int>(this->Vertices.size());
    }

    //
    // NumEdges
    //
    // Returns the # of edges currently in the graph.
    //
    int NumEdges() const {
        int count = 0;
        //
        // loop through the adjacency matrix and count how many
        // edges currently exist:
        //
        for (auto i : adjList) {
            count = count + i.second.size();
        }
        return count;
    }

    //
    // addVertex
    //
    // Adds the vertex v to the graph if there's room, and if so
    // returns true.  If the graph is full, or the vertex already
    // exists in the graph, then false is returned.
    //
    bool addVertex(VertexT v) {
        //
        // is the vertex already in the graph?  If so, we do not
        // insert again otherwise Vertices may fill with duplicates:
        //
        if (Vertices.count(v) == 1) {
            return false;
        }
        this->vertexs.push_back(v);
        //
        // if we get here, vertex does not exist so insert.  Where
        // we insert becomes the rows and col position for this
        // vertex in the adjacency matrix.
        //
        this->Vertices.insert(v);
        return true;
    }

    //
    // addEdge
    //
    // Adds the edge (from, to, weight) to the graph, and returns
    // true.  If the vertices do not exist or for some reason the
    // graph is full, false is returned.
    //
    // NOTE: if the edge already exists, the existing edge weight
    // is overwritten with the new edge weight.
    //
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        //
        // we need to search the Vertices and find the position
        // of each vertex; this will denote the row and col to
        // access in the adjacency matrix:
        if (Vertices.count(from) == 0) {  // not find
            return false;
        }
        if (Vertices.count(to) == 0) {   // not find from set
            return false;
        }
        // fing the edge
        adjList[from][to] = weight;
        return true;
    }

    //
    // getWeight
    //
    // Returns the weight associated with a given edge.  If
    // the edge exists, the weight is returned via the reference
    // parameter and true is returned.  If the edge does not
    // exist, the weight parameter is unchanged and false is
    // returned.
    //
    bool getWeight(VertexT from, VertexT to, WeightT& weight) {
        // find if we have the vertex
        if (Vertices.count(from) == 0) {
            return false;
        }
        if (Vertices.count(to) == 0) {
            return false;
        }
        // find the edge and return true with update the weight
        if (adjList[from].count(to) == 1) {
            weight = adjList[from][to];
            return true;
        }
        return false;
    }
    //
    // neighbors
    //
    // Returns a set containing the neighbors of v, i.e. all
    // vertices that can be reached from v along one edge.
    // Since a set is returned, the neighbors are returned in
    // sorted order; use foreach to iterate through the set.
    //
    set<VertexT> neighbors(VertexT v) {
        set<VertexT> S;
        // find if we have the vertex
        if (Vertices.count(v) == 0) {
            return S;
        }
        // insert the neighbors into set
        for (auto i : adjList[v]) {
            S.insert(i.first);
        }
        return S;
    }
    //
    // getVertices
    //
    // Returns a vector containing all the vertices currently in
    // the graph.
    //
    vector<VertexT> getVertices() const {
        return this->vertexs;  // returns a copy:
    }


    // this is the equal operator, it copy one graph class
    // into another
    graph& operator=(const graph& other) {
        if (this == &other) {
            return *this;
        }
        // copy three private variables
        adjList = other.adjList;
        Vertices = other.Vertices;
        vertexs = other.vertexs;
        return *this;
    }

    // this is the copy constructor, used for functions parameter;
    graph(const graph& other) {
        adjList.clear();
        Vertices.clear();
        vertexs.clear();
        // copy three private variables
        adjList = other.adjList;
        Vertices = other.Vertices;
        vertexs = other.vertexs;
    }
    //
    // dump
    //
    // Dumps the internal state of the graph for debugging purposes.
    //
    // Example:
    //    graph<string,int>  G(26);
    //    ...
    //    G.dump(cout);  // dump to console
    //
    void dump(ostream& output) const {
        output << "***************************************************" << endl;
        output << "********************* GRAPH ***********************" << endl;
        output << "**Num vertices: " << this->NumVertices() << endl;
        output << "**Num edges: " << this->NumEdges() << endl;
        output << endl;
        output << "**Vertices:" << endl;
        int a = 0;
        for (auto i : Vertices) {
            output << " " << a << ". " << i << endl;
            ++a;
        }
        output << endl;
        output << "**Edges:" << endl;
        for (auto i : adjList) {
            output << i.first << ": ";
            for (auto j : i.second) {
                output << "(" << i.first << ",";
                output << j.first << "," << j.second << ") ";
            }
            output << endl;
        }
        output << "**************************************************" << endl;
    }
};
   
