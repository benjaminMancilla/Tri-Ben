#ifndef PREDICATES_HPP
#define PREDICATES_HPP

#include "triangle.hpp"
#include <limits>
#include <cmath>

template <typename T> class Predicates {
    public:
    T epsilon = std::numeric_limits<T>::epsilon() * 1e8;
    bool in_circle(const Triangle<T>* t, const Vertex<T>* v) const;
    int in_triangle(const Triangle<T>* t, const Vertex<T>* v);
    bool orientation(const Triangle<T>* t) const;
    T orientationPoint(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* P) const;
    bool is_between(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* P);
    bool is_between_loose(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* P);
    std::vector<Vertex<T>*> find_collinear_points(const Vertex<T>* A, const Vertex<T>* B, const std::vector<Vertex<T>*>* points);
    bool doIntersect(Vertex<T>* p1, Vertex<T>* q1, Vertex<T>* p2, Vertex<T>* q2, Vertex<T>* intersection = nullptr);
    bool doIntersect_loose(Vertex<T>* p1, Vertex<T>* q1, Vertex<T>* p2, Vertex<T>* q2, Vertex<T>* intersection = nullptr);
    //bool is_left_of(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* P) const;
    

};

template <typename T> bool Predicates<T>::in_circle(const Triangle<T>* t, const Vertex<T>* v) const {
    //Calcula el determinante de la matriz de la circunferencia
    //de un triangulo y un punto
    // |x1 y1 x1^2 + y1^2 1|
    // |x2 y2 x2^2 + y2^2 1|
    // |x3 y3 x3^2 + y3^2 1|
    // |x4 y4 x4^2 + y4^2 1|
    // Si el determinante es mayor a 0, el punto esta dentro de la circunferencia
    // del triangulo
    /*
    function inCircle (ax, ay, bx, by, cx, cy, dx, dy) {
    let ax_ = ax-dx;
    let ay_ = ay-dy;
    let bx_ = bx-dx;
    let by_ = by-dy;
    let cx_ = cx-dx;
    let cy_ = cy-dy;
    return (
        (ax_*ax_ + ay_*ay_) * (bx_*cy_-cx_*by_) -
        (bx_*bx_ + by_*by_) * (ax_*cy_-cx_*ay_) +
        (cx_*cx_ + cy_*cy_) * (ax_*by_-bx_*ay_)
    ) > 0;
}
    */
    T ax = t->vertices[0]->x_coord;
    T ay = t->vertices[0]->y_coord;
    T bx = t->vertices[1]->x_coord;
    T by = t->vertices[1]->y_coord;
    T cx = t->vertices[2]->x_coord;
    T cy = t->vertices[2]->y_coord;
    T dx = v->x_coord;
    T dy = v->y_coord;
    T ax_ = ax - dx;
    T ay_ = ay - dy;
    T bx_ = bx - dx;
    T by_ = by - dy;
    T cx_ = cx - dx;
    T cy_ = cy - dy;

    T determinant = (ax_ * ax_ + ay_ * ay_) * (bx_ * cy_ - cx_ * by_) -
                    (bx_ * bx_ + by_ * by_) * (ax_ * cy_ - cx_ * ay_) +
                    (cx_ * cx_ + cy_ * cy_) * (ax_ * by_ - bx_ * ay_);

    bool orientation = this->orientation(t);
    if (!orientation) {
        return determinant < -1*epsilon;
    }
    return determinant > epsilon;
}

template <typename T> int Predicates<T>::in_triangle(const Triangle<T>* t, const Vertex<T>* v) {
    //Calcula el determinante de la matriz de un triangulo y un punto
    // |x1 y1 1|
    // |x2 y2 1|
    // |x3 y3 1|
    // |x4 y4 1|
    // Si el determinante es mayor a 0, el punto esta dentro del triangulo

    T D1, D2, D3;
    bool has_neg, has_pos;


    D1 = orientationPoint(v, t->vertices[0], t->vertices[1]);
    D2 = orientationPoint(v, t->vertices[1], t->vertices[2]);
    D3 = orientationPoint(v, t->vertices[2], t->vertices[0]);

    auto distance = [&](Vertex<T>* A, Vertex<T>* B) -> T {
        return std::sqrt(std::pow(A->x_coord - B->x_coord, 2) + std::pow(A->y_coord - B->y_coord, 2));
    };
    T distance1 = distance(t->vertices[0], t->vertices[1]);
    T distance2 = distance(t->vertices[1], t->vertices[2]);
    T distance3 = distance(t->vertices[2], t->vertices[0]);

    auto is_almost_zero = [&](T value, T scale_factor) -> bool {
        //T result = epsilon; //* pow(scale_factor, 3.5);//se queda asi mientras

        return std::abs(value) < epsilon;
    };
    
    has_neg = (D1 < -epsilon) || (D2 < -epsilon) || (D3 < -epsilon);
    has_pos = (D1 > epsilon) || (D2 > epsilon) || (D3 > epsilon);

    bool inside_condition = !(has_neg && has_pos);

    if (inside_condition) {
        if (is_almost_zero(D1, distance1)) {
            //cout << "D1: " << D1 << endl;
            return 2; 
        }
        if(is_almost_zero(D2, distance2)){
            //cout << "D2: " << D2 << endl;
            return 3;
        }
        if(is_almost_zero(D3, distance3)){
            //cout << "D3: " << D3 << endl;
            return 4;
        }
        return 1;
    }

    return 0;
}

template <typename T> bool Predicates<T>::orientation(const Triangle<T>* t) const {
    //Calcula el determinante de la matriz de un triangulo
    // |x1 y1 1|
    // |x2 y2 1|
    // |x3 y3 1|
    // Si el determinante es mayor a 0, el triangulo es orientado en sentido horario
    T Ax = t->vertices[0]->x_coord;
    T Ay = t->vertices[0]->y_coord;
    T Bx = t->vertices[1]->x_coord;
    T By = t->vertices[1]->y_coord;
    T Cx = t->vertices[2]->x_coord;
    T Cy = t->vertices[2]->y_coord;

    T det = ((Ax - Cx) * (By - Cy) - (Ay - Cy) * (Bx - Cx));
    return bool(det > 0);
}

//template <typename T> bool Predicates<T>::is_left_of(const Vertex<T> &A, const Vertex<T> &B, const Vertex<T> &P) const {
//    return ((B->x_coord - A->x_coord) * (P->y_coord - A->y_coord) - (B->y_coord - A->y_coord) * (P->x_coord - A->x_coord)) > 0;
//}

template <typename T> T Predicates<T>::orientationPoint(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* C) const {
    return ((A->x_coord - C->x_coord) * (B->y_coord - C->y_coord) - (A->y_coord - C->y_coord) * (B->x_coord - C->x_coord));
    

    
}

const double epsilon = 1e-9;  // Puedes ajustar el valor según la precisión que necesites

template <typename T>
bool Predicates<T>::is_between(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* P) {
    if (std::abs(A->x_coord - B->x_coord) < epsilon) {  // Segmento vertical
        return std::abs(P->x_coord - A->x_coord) < epsilon &&
               std::min(A->y_coord, B->y_coord) - epsilon <= P->y_coord &&
               P->y_coord <= std::max(A->y_coord, B->y_coord) + epsilon;
    } else if (std::abs(A->y_coord - B->y_coord) < epsilon) {  // Segmento horizontal
        return std::abs(P->y_coord - A->y_coord) < epsilon &&
               std::min(A->x_coord, B->x_coord) - epsilon <= P->x_coord &&
               P->x_coord <= std::max(A->x_coord, B->x_coord) + epsilon;
    } else {  // Caso general
        return std::min(A->x_coord, B->x_coord) - epsilon <= P->x_coord &&
               P->x_coord <= std::max(A->x_coord, B->x_coord) + epsilon &&
               std::min(A->y_coord, B->y_coord) - epsilon <= P->y_coord &&
               P->y_coord <= std::max(A->y_coord, B->y_coord) + epsilon;
    }
}

template <typename T>
bool Predicates<T>::is_between_loose(const Vertex<T>* A, const Vertex<T>* B, const Vertex<T>* P) {
    if (std::abs(A->x_coord - B->x_coord) < epsilon) {  // Segmento vertical
        return std::abs(P->x_coord - A->x_coord) < epsilon &&
               std::min(A->y_coord, B->y_coord) + epsilon < P->y_coord &&
               P->y_coord < std::max(A->y_coord, B->y_coord) - epsilon;
    } else if (std::abs(A->y_coord - B->y_coord) < epsilon) {  // Segmento horizontal
        return std::abs(P->y_coord - A->y_coord) < epsilon &&
               std::min(A->x_coord, B->x_coord) + epsilon < P->x_coord &&
               P->x_coord < std::max(A->x_coord, B->x_coord) - epsilon;
    } else {  // Caso general
        return std::min(A->x_coord, B->x_coord) + epsilon < P->x_coord &&
               P->x_coord < std::max(A->x_coord, B->x_coord) - epsilon &&
               std::min(A->y_coord, B->y_coord) + epsilon < P->y_coord &&
               P->y_coord < std::max(A->y_coord, B->y_coord) - epsilon;
    }
}




template <typename T> 
std::vector<Vertex<T>*> Predicates<T>::find_collinear_points(const Vertex<T>* A, const Vertex<T>* B, const std::vector<Vertex<T>*>* points) {
    std::vector<Vertex<T>*> result;

    for (const auto p : *points) {
        //cout << "p: " << p->x_coord << " - " << p->y_coord << endl;
        if (*p == *A || *p == *B){
            //cout << "p is A or B" << endl;
            result.push_back(p);
            continue;
        }
        if (p->x_coord < std::min(A->x_coord, B->x_coord) || p->x_coord > std::max(A->x_coord, B->x_coord)){
            //cout << "x out of range" << endl;
            continue;
        };
        if (p->y_coord < std::min(A->y_coord, B->y_coord) || p->y_coord > std::max(A->y_coord, B->y_coord)){
            //cout << "y out of range" << endl;
            continue;
        };

        if ((orientationPoint(A, B, p) == 0) && is_between_loose(A, B, p)) {
            //cout << "p is collinear" << endl;
            result.push_back(p);
            continue;
        }
        //cout << "p is not collinear" << endl;
    }
    return result;
}

template <typename T> bool Predicates<T>::doIntersect(Vertex<T>* p1, Vertex<T>* q1, Vertex<T>* p2, Vertex<T>* q2, Vertex<T>* intersection)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientationPoint(p1, q1, p2);
    int o2 = orientationPoint(p1, q1, q2);
    int o3 = orientationPoint(p2, q2, p1);
    int o4 = orientationPoint(p2, q2, q1);

    // General case: lines intersect if they have different
    // orientations
    if (o1 != o2 && o3 != o4) {
        // Compute intersection point
        T a1 = q1->y_coord - p1->y_coord;
        T b1 = p1->x_coord - q1->x_coord;
        T c1 = a1 * p1->x_coord + b1 * p1->y_coord;

        T a2 = q2->y_coord - p2->y_coord;
        T b2 = p2->x_coord - q2->x_coord;
        T c2 = a2 * p2->x_coord + b2 * p2->y_coord;

        double determinant = a1 * b2 - a2 * b1;
        Vertex<T> intersection_aux;
        //cout << "Intersecting : " << *p1 << " - " << *q1 << " with " << *p2 << " - " << *q2 << endl;

        if (determinant != 0) {
            intersection_aux.x_coord
                = (c1 * b2 - c2 * b1) / determinant;
            intersection_aux.y_coord
                = (a1 * c2 - a2 * c1) / determinant;
            //cout << "Intersection: " << intersection_aux << endl;
            if(intersection != nullptr){
                intersection->x_coord = intersection_aux.x_coord;
                intersection->y_coord = intersection_aux.y_coord;
            }
            if(is_between_loose(p1, q1, &intersection_aux) && is_between(p2, q2, &intersection_aux)){
                cout << "Intersect" << endl;
                return true;
            }
            else{
                cout << "Intersect but not between" << endl;
                cout << "P1: " << *p1 << endl;
                cout << "Q1: " << *q1 << endl;
                cout << "P2: " << *p2 << endl;
                cout << "Q2: " << *q2 << endl; 
                cout << "Intersection: " << intersection_aux << endl;


                return false;
            }
        }
    }
    /*
    cout << "Special case intersection" << endl;
    cout << "P1: " << *p1 << endl;
    cout << "Q1: " << *q1 << endl;
    cout << "P2: " << *p2 << endl;
    cout << "Q2: " << *q2 << endl;
    */

    // Special Cases: check if the lines are collinear and
    // overlap
    if (o1 == 0 && is_between(p1, p2, q1))
        return true;
    if (o2 == 0 && is_between(p1, q2, q1))
        return true;
    // Lines do not intersect in any case
    return false;
}

template <typename T> bool Predicates<T>::doIntersect_loose(Vertex<T>* p1, Vertex<T>* q1, Vertex<T>* p2, Vertex<T>* q2, Vertex<T>* intersection)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientationPoint(p1, q1, p2);
    int o2 = orientationPoint(p1, q1, q2);
    int o3 = orientationPoint(p2, q2, p1);
    int o4 = orientationPoint(p2, q2, q1);

    // General case: lines intersect if they have different
    // orientations
    if (o1 != o2 && o3 != o4) {
        // Compute intersection point
        T a1 = q1->y_coord - p1->y_coord;
        T b1 = p1->x_coord - q1->x_coord;
        T c1 = a1 * p1->x_coord + b1 * p1->y_coord;

        T a2 = q2->y_coord - p2->y_coord;
        T b2 = p2->x_coord - q2->x_coord;
        T c2 = a2 * p2->x_coord + b2 * p2->y_coord;

        double determinant = a1 * b2 - a2 * b1;
        Vertex<T> intersection_aux;
        //cout << "Intersecting : " << *p1 << " - " << *q1 << " with " << *p2 << " - " << *q2 << endl;
        //cout << "Determinant: " << determinant << endl;
        if (determinant != 0) {
            intersection_aux.x_coord
                = (c1 * b2 - c2 * b1) / determinant;
            intersection_aux.y_coord
                = (a1 * c2 - a2 * c1) / determinant;
            //cout << "Intersection: " << intersection_aux << endl;
            if(intersection != nullptr){
                intersection->x_coord = intersection_aux.x_coord;
                intersection->y_coord = intersection_aux.y_coord;
            }
            if(is_between_loose(p1, q1, &intersection_aux) && is_between_loose(p2, q2, &intersection_aux)){//Deje OR pero en verdad hay que dejar el isBetween mas flexible
                //cout << "Intersect" << endl;
                return true;
            }
            else{
                //cout << "Intersect but not between" << endl;
                return false;
            }
        }
    }

    // Special Cases: check if the lines are collinear and
    // overlap
    if (o1 == 0 && is_between_loose(p1, p2, q1))
        return true;
    if (o2 == 0 && is_between_loose(p1, q2, q1))
        return true;
    // Lines do not intersect in any case
    return false;
}



#endif // PREDICATES_HPP