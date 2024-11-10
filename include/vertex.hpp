#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <iostream>

using namespace std;

template <typename T> class Vertex {
public:
    //Tuple of coordinates
    T x_coord;
    T y_coord;
    int id = 0;
    Vertex(T x, T y);
    Vertex();
    ~Vertex();
    Vertex(const Vertex<T> &v);
    bool operator==(const Vertex<T> &v) const;
    bool operator!=(const Vertex<T> &v) const;
};

template <typename T> Vertex<T>::Vertex(T x, T y) {
    x_coord = x;
    y_coord = y;
}

template <typename T> Vertex<T>::Vertex() {
    x_coord = 0;
    y_coord = 0;
}

template <typename T> Vertex<T>::~Vertex() = default;

template <typename T> Vertex<T>::Vertex(const Vertex<T> &v) {
    x_coord = v.x_coord;
    y_coord = v.y_coord;
}

//Usar predicados para comparar
template <typename T> bool Vertex<T>::operator==(const Vertex<T> &v) const {
    return (x_coord == v.x_coord && y_coord == v.y_coord);
}

template <typename T> bool Vertex<T>::operator!=(const Vertex<T> &v) const {
    return (x_coord != v.x_coord || y_coord != v.y_coord);
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const Vertex<T> &v) {
    os << "(" << v.x_coord << ", " << v.y_coord << ")";
    return os;
}



#endif // VERTEX_HPP
