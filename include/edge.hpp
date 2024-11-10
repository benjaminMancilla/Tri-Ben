#ifndef EDGE_HPP
#define EDGE_HPP

#include <vector>
#include "triangle.hpp"

template <typename T> 
class Edge {
public:
    using TrianglePair = std::pair<Triangle<T>*, int>;
    std::vector<TrianglePair> adjacent_triangles;
    std::vector<Vertex<T>*> vertices;
    bool constrained = false;
    bool special = false;
    int original_triangle = 0;

    Edge();
    ~Edge();
    Edge(const Edge<T> &e);
    Edge(const Vertex<T>* v1, const Vertex<T>* v2);
    void add_adjacent_triangle(TrianglePair t);
    void replace_adjacent_triangle(const Triangle<T> *old, const TrianglePair& new_t);
    Triangle<T>* getLeftTriangle() const;
    Triangle<T>* getRightTriangle() const;
    int getLeftTriangleEdgeId() const;
    int getRightTriangleEdgeId() const;
    Triangle<T>* getOtherTriangle(const Triangle<T>* t) const;
    bool isBoundaryEdge() const;
    int get_index_from_triangle(const Triangle<T>* t) const;
    void set_vertices(Vertex<T>* v1, Vertex<T>* v2);
    T length() const;
};

template <typename T> 
Triangle<T>* Edge<T>::getLeftTriangle() const {
    return adjacent_triangles[original_triangle].first;
}

template <typename T>
int Edge<T>::get_index_from_triangle(const Triangle<T>* t) const {
    if (adjacent_triangles[0].first == t) {
        return adjacent_triangles[0].second;
    } else {
        return adjacent_triangles[1].second;
    }
}

template <typename T> 
Triangle<T>* Edge<T>::getRightTriangle() const {
    return adjacent_triangles[!original_triangle].first;
}

template <typename T> 
int Edge<T>::getLeftTriangleEdgeId() const {
    return adjacent_triangles[original_triangle].second;
}

template <typename T> 
int Edge<T>::getRightTriangleEdgeId() const {
    return adjacent_triangles[!original_triangle].second;
}

template <typename T> 
Edge<T>::Edge()
    :adjacent_triangles(2)
{

}

template <typename T>
Edge<T>::Edge(const Vertex<T>* v1, const Vertex<T>* v2) 
    : adjacent_triangles(2) {
        vertices.push_back(const_cast<Vertex<T>*>(v1));
        vertices.push_back(const_cast<Vertex<T>*>(v2));
        constrained = true;
}


template <typename T> 
Edge<T>::~Edge() {
    adjacent_triangles.clear();
}

template <typename T> 
Edge<T>::Edge(const Edge<T> &e) {
    adjacent_triangles = e.adjacent_triangles;
    original_triangle = e.original_triangle;
    vertices = e.vertices;
    constrained = e.constrained;
    original_triangle = e.original_triangle;
}

template <typename T> 
void Edge<T>::add_adjacent_triangle(TrianglePair t) {
    if (adjacent_triangles[0].first == nullptr) {
        adjacent_triangles[0] = t;
    } else {
        adjacent_triangles[1] = t;
    }
}

template <typename T> 
void Edge<T>::replace_adjacent_triangle(const Triangle<T> *old, const TrianglePair& new_t) {
    if (adjacent_triangles[0].first == old) {
        adjacent_triangles[0] = new_t;
    } else {
        adjacent_triangles[1] = new_t;
    }
}

template <typename T> 
bool Edge<T>::isBoundaryEdge() const {
    return adjacent_triangles.size() == 1;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const Edge<T> &e) {
    if(e.adjacent_triangles[0].first == nullptr){
        os << "Constrained Edge: " << e.vertices[0] << " - " << e.vertices[1];
        return os;
    }
    Triangle<T>* t1 = e.adjacent_triangles[0].first;
    int e1 = e.adjacent_triangles[0].second;
    Vertex<T>* v1 = t1->vertices[e1];
    Vertex<T>* v2 = t1->vertices[(e1 + 1) % 3];
    os << *v1 << " - " << *v2;
    return os;
}

template <typename T>
Triangle<T>* Edge<T>::getOtherTriangle(const Triangle<T>* t) const {
    if (adjacent_triangles[0].first == t) {
        return adjacent_triangles[1].first;
    } else {
        return adjacent_triangles[0].first;
    }
}

template <typename T>
void Edge<T>::set_vertices(Vertex<T>* v1, Vertex<T>* v2) {
    vertices.push_back(v1);
    vertices.push_back(v2);
}

template <typename T>
T Edge<T>::length() const {
    return sqrt(pow(vertices[0]->x_coord - vertices[1]->x_coord, 2) + pow(vertices[0]->y_coord - vertices[1]->y_coord, 2));
}



/*
void Edge::determineAdjacentTriangles(int apexVertexLeftTriangleId)
{
    this->correctOrientation =
            adjacentTrianglesInfo[0].first->vertices[minus1mod3[adjacentTrianglesInfo[0].second]]->id ==
            apexVertexLeftTriangleId;
}
*/

#endif // EDGE_HPP

