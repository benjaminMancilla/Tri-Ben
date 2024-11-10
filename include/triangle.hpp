#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include "vertex.hpp"
#include <vector>

template <typename T> class Edge;

template <typename T> class Triangle {
public:
    std::vector<Vertex<T>*> vertices;
    std::vector<Edge<T>*> edges;
    bool boundTriangle = false;
    bool in_da_wae = true;
    int debug_number;

    Triangle(Vertex<T> *v1, Vertex<T> *v2, Vertex<T> *v3);
    Triangle(std::vector<Vertex<T>*> vertices);
    Triangle();
    ~Triangle();
    Triangle(const Triangle<T> &t);
    bool operator==(const Triangle<T> &t) const;
    bool operator!=(const Triangle<T> &t) const;
    void set_edges(Edge<T> *e1, Edge<T> *e2, Edge<T> *e3);
    bool contains_vertex(const Vertex<T>* v) const;
    bool isBoundaryTriangle() const;

};

template <typename T> Triangle<T>::Triangle(Vertex<T> *v1, Vertex<T> *v2, Vertex<T> *v3) {
    vertices.push_back(v1);
    vertices.push_back(v2);
    vertices.push_back(v3);
}

template <typename T> Triangle<T>::Triangle(std::vector<Vertex<T>*> vertices) {
    this->vertices = vertices;
}

template <typename T> Triangle<T>::Triangle() {
    vertices = std::vector<Vertex<T>*>(3);
    edges = std::vector<Edge<T>*>(3);
}


template <typename T>
Triangle<T>::~Triangle() {
    vertices.clear();
    edges.clear();
};
//{
 //   for (Edge<T>* edge : edges) {
  //      if (edge != nullptr) {
   //         delete edge;
    //    }
    //}
//}

template <typename T> Triangle<T>::Triangle(const Triangle<T> &t) {
    vertices = t.vertices;
    edges = t.edges;
}

template <typename T> bool Triangle<T>::operator==(const Triangle<T> &t) const {
    //Ojito, no deberian haber triangulos con los mismos vertices en distinto
    //orden, pero si ocurre algun bucle puede ser aca
    for (int i = 0; i < 3; i++) {
        if (vertices[i] != t.vertices[i]) {
            return false;
        }
    }
}

template <typename T> bool Triangle<T>::operator!=(const Triangle<T> &t) const {
    return !(*this == t);
}

template <typename T> void Triangle<T>::set_edges(Edge<T> *e1, Edge<T> *e2, Edge<T> *e3) {
    edges.push_back(e1);
    edges.push_back(e2);
    edges.push_back(e3);
}

template <typename T> bool Triangle<T>::contains_vertex(const Vertex<T>* v) const {
    for (int i = 0; i < 3; i++) {
        if (*(vertices[i]) == *v) {
            return true;
        }
    }
    return false;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const Triangle<T> &t) {
    Vertex<T>* v1_p = t.vertices[0];
    Vertex<T>* v2_p = t.vertices[1];
    Vertex<T>* v3_p = t.vertices[2];
    Vertex<T> v1 = *v1_p;
    Vertex<T> v2 = *v2_p;
    Vertex<T> v3 = *v3_p;
    os << v1 << " - " << v2 << " - " << v3;
    return os;
}

template <typename T> bool Triangle<T>::isBoundaryTriangle() const {
    for (auto vertex : vertices) {
        if (vertex->id < 0) {
            return true;
        }
    }
    return false;
}



#endif // TRIANGLE_HPP