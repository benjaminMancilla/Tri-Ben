#ifndef DELAUNAY_TRIANGULATION_HPP
#define DELAUNAY_TRIANGULATION_HPP

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include "edge.hpp"
#include "predicates.hpp"


template <typename T> class DelaunayTriangulation {
public:
    using TrianglePair = std::pair<Triangle<T>*, int>;
    std::vector<Edge<T>*> mesh_edges;
    std::vector<Triangle<T>*> mesh_triangles;
    std::vector<Vertex<T>*> mesh_vertices;
    std::vector<Edge<T>*> mesh_edges_debug;
    std::vector<Triangle<T>*> mesh_triangles_debug;
    std::vector<Vertex<T>*> mesh_vertices_debug;
    std::vector<Edge<T>*> constraid_edges;
    int n_cache = 0;

    DelaunayTriangulation();
    ~DelaunayTriangulation();
    void initVertex_and_ConsEdges(std::vector<Vertex<T>> &vertices, std::vector<Edge<T>> &edges);
    void boundingTriangle();
    void insertVertex(Vertex<T> &v);
    void legalizeEdge(Triangle<T>* original_triangle, Edge<T>* illegal_edge, Vertex<T>* v);
    void triangulate(int n = -1);
    Triangle<T>* locate(Vertex<T>* Pr, int* insert_type=nullptr);
    void print_input_vertex();
    void print_edges();
    void print_triangles();
    std::vector<Triangle<T>> getTriangles();
    std::vector<Edge<T>> getEdges();
    std::vector<Vertex<T>> getVertices();
    void find_in_path_triangles(std::vector<std::pair<Triangle<T>*, Edge<T>*>>* triangles_to_flip, Vertex<T>* Pi, Vertex<T>* Pj);
    void this_is_the_way(std::vector<std::pair<Triangle<T>*, Edge<T>*>>* triangles_to_flip, Vertex<T>* Pi, Vertex<T>* Pj);
    bool da_way_flip(Triangle<T>* PiPjPr, Edge<T>* PiPj, 
        std::pair<std::vector<std::pair<Triangle<T>*, Edge<T>*>>*, int> triangles_to_update = {{}, 0}, Vertex<T>* Pi = nullptr, Vertex<T>* Pj = nullptr);
    void insert_constrained_edges();
    void remove_super_triangle();
    void add_constrainedEdge(Vertex<T>* v1, Vertex<T>* v2);
    std::pair<std::vector<Triangle<T>*>,Edge<T>*> LEPP_search(Triangle<T>* t);
};

template <typename T> DelaunayTriangulation<T>::DelaunayTriangulation() {
    mesh_edges = std::vector<Edge<T>*>();
    mesh_triangles = std::vector<Triangle<T>*>();
    mesh_vertices = std::vector<Vertex<T>*>();
    constraid_edges = std::vector<Edge<T>*>();
}

template <typename T> DelaunayTriangulation<T>::~DelaunayTriangulation() {
    mesh_edges.clear();
    mesh_triangles.clear();
    mesh_vertices.clear();
    constraid_edges.clear();
}

template <typename T>
void DelaunayTriangulation<T>::initVertex_and_ConsEdges(std::vector<Vertex<T>>& vertices, std::vector<Edge<T>>& edges) {
    std::sort(vertices.begin(), vertices.end(), [](Vertex<T>& v1, Vertex<T>& v2) {
        return v1.x_coord < v2.x_coord || (v1.x_coord == v2.x_coord && v1.y_coord < v2.y_coord);
    });


    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());

    std::random_shuffle(vertices.begin(), vertices.end());

    mesh_vertices.clear();
    mesh_triangles.clear();
    mesh_edges.clear();
    constraid_edges.clear();

    for (size_t i = 0; i < vertices.size(); i++) {
        auto vertex = new Vertex<T>(vertices[i]);
        vertex->id = static_cast<int>(i);
        mesh_vertices.push_back(vertex);
    }

    vertices.clear();

    for (auto& edge : edges) {
        Vertex<T>* v1 = nullptr;
        Vertex<T>* v2 = nullptr;

        for (auto& mesh_vertex : mesh_vertices) {
            if (mesh_vertex->x_coord == edge.vertices[0]->x_coord && mesh_vertex->y_coord == edge.vertices[0]->y_coord) {
                v1 = mesh_vertex;
            }
            if (mesh_vertex->x_coord == edge.vertices[1]->x_coord && mesh_vertex->y_coord == edge.vertices[1]->y_coord) {
                v2 = mesh_vertex;
            }
        }

        if (v1 != nullptr && v2 != nullptr) {
            auto restricted_edge = new Edge<T>(v1, v2);
            constraid_edges.push_back(restricted_edge);
        } else {
            std::cerr << "Error: No se encontraron los vértices para el edge restringido en la malla." << std::endl;
        }
    }

    edges.clear();
}

template<typename T> void DelaunayTriangulation<T>::add_constrainedEdge(Vertex<T>* v1, Vertex<T>* v2) {
    auto edge = new Edge<T>(v1, v2);
    constraid_edges.push_back(edge);
}


template <typename T> void DelaunayTriangulation<T>::boundingTriangle() {
    T min_x = mesh_vertices[0]->x_coord;
    T min_y = mesh_vertices[0]->y_coord;
    T max_x = min_x;
    T max_y = min_y;

    for (auto &v : mesh_vertices) {
        if (v->x_coord < min_x) {
            min_x = v->x_coord;
        }
        if (v->y_coord < min_y) {
            min_y = v->y_coord;
        }
        if (v->x_coord > max_x) {
            max_x = v->x_coord;
        }
        if (v->y_coord > max_y) {
            max_y = v->y_coord;
        }
    }

    T dx = max_x - min_x;
    T dy = max_y - min_y;
    T mid_x = (max_x + min_x) / 2;
    T mid_y = (max_y + min_y) / 2;
    T delta = std::max(dx, dy);
    double multiplier = 100;

    auto Pi = new Vertex<T>({mid_x - multiplier * delta, mid_y - multiplier*delta/sqrt(3)});
    Pi->id = -3;
    auto Pj = new Vertex<T>({mid_x, mid_y + multiplier * delta});
    Pj->id = -2;
    auto Pk = new Vertex<T>({mid_x + multiplier * delta, mid_y - multiplier*delta/sqrt(3)});
    Pk->id = -1;
    mesh_vertices.push_back(Pi);
    mesh_vertices.push_back(Pj);
    mesh_vertices.push_back(Pk);

    auto PiPj = new Edge<T>();
    PiPj->set_vertices(Pi, Pj);
    auto PjPk = new Edge<T>();
    PjPk->set_vertices(Pj, Pk);
    auto PkPi = new Edge<T>();
    PkPi->set_vertices(Pk, Pi);

    auto PiPjPk = new Triangle<T>(Pi, Pj, Pk);
    PiPjPk->set_edges(PiPj, PjPk, PkPi);
    PiPjPk->debug_number = -1;

    PiPj->add_adjacent_triangle(TrianglePair(PiPjPk, 0));
    PjPk->add_adjacent_triangle(TrianglePair(PiPjPk, 1));
    PkPi->add_adjacent_triangle(TrianglePair(PiPjPk, 2));

    mesh_edges.push_back(PiPj);
    mesh_edges.push_back(PjPk);
    mesh_edges.push_back(PkPi);

    mesh_triangles.push_back(PiPjPk);


}

template <typename T> std::pair<std::vector<Triangle<T>*>,Edge<T>*> DelaunayTriangulation<T>::LEPP_search(Triangle<T>* t_0) {
    //Initialize LEPP path and terminal edge
    std::vector<Triangle<T>*> LEPP_path;
    Edge<T>* terminal_edge = nullptr;
    std::unordered_set<Triangle<T>*> visited_triangles;

    if (t_0 == nullptr) {
        std::cerr << "Error: Initial triangle is nullptr" << std::endl;
        return {LEPP_path, nullptr};
    }

    //Add t_0 to LEPP path
    LEPP_path.push_back(t_0);
    visited_triangles.insert(t_0);
    //We add triangles to the path until we find a terminal edge
    while(terminal_edge == nullptr){
        //Get the current triangle and its edges
        Triangle<T>* currentT = LEPP_path.back();

        //Get the edge with the longest length
        Edge<T>* longest_edge = nullptr;
        T longest_length = 0;
        for (int i = 0; i < 3; i++) {
            auto edge = currentT->edges[i];
            if (edge != nullptr && edge->length() > longest_length) {
                longest_edge = edge;
                longest_length = edge->length();
            }
        }

        if (longest_edge == nullptr) {
            std::cerr << "Error: Longest edge is nullptr" << std::endl;
            return {LEPP_path, nullptr};
        }

        //Check if the edge is boundary
        if(longest_edge->isBoundaryEdge()){
            terminal_edge = longest_edge;
            break;
        }

        //Get the other triangle of the longest edge
        auto next_triangle = longest_edge->getOtherTriangle(currentT);
        if (next_triangle == nullptr) {
            std::cerr << "Error: No neighboring triangle found" << std::endl;
            terminal_edge = longest_edge;
            break;
        }

        if(next_triangle->isBoundaryTriangle()){
            terminal_edge = longest_edge;
            break;
        }


        //Check if the next triangle is already in the path
        if (visited_triangles.count(next_triangle) > 0) {
            terminal_edge = longest_edge;
            break;
        }

        //Add the next triangle to the path and continue
        LEPP_path.push_back(next_triangle);
        visited_triangles.insert(next_triangle);
        
    }

    return {LEPP_path, terminal_edge};
    
}


template <typename T> Triangle<T>* DelaunayTriangulation<T>::locate(Vertex<T>* Pr, int* insert_type) {
    Predicates<T> pred;
    Triangle<T>* current_triangle = mesh_triangles[0];
    while (true) {
        int in_triangle = pred.in_triangle(current_triangle, Pr);
        if (in_triangle == 1) {
            if (insert_type != nullptr){
                *insert_type = -1;
            }
            return current_triangle;

        } else if (in_triangle == 2 || in_triangle == 3 || in_triangle == 4) {
            if(insert_type != nullptr){
                *insert_type = in_triangle-2;
            }
            return current_triangle;

        } else {
            Edge<T>* edge = nullptr;
            //True is CounterClockwise
            //False is Clockwise
            bool triangle_orientation = pred.orientation(current_triangle);
            bool desired_orientation = !triangle_orientation;
            Triangle<T>* next_triangle;
            for (int i = 0; i < 3; i++) {
                T edge_pr_orientation = pred.orientationPoint(Pr, current_triangle->vertices[i], current_triangle->vertices[(i + 1) % 3]);
                if (edge_pr_orientation > 0 && desired_orientation) {
                    edge = current_triangle->edges[i];
                    next_triangle = edge->getOtherTriangle(current_triangle);
                    if (next_triangle == nullptr) {
                        cout << "Error: No triangle found in locate" << endl;
                        return nullptr;
                    }
                    break;
                }
                if (edge_pr_orientation < 0 && !desired_orientation) {
                    edge = current_triangle->edges[i];
                    edge = current_triangle->edges[i];
                    next_triangle = edge->getOtherTriangle(current_triangle);
                    if (next_triangle == nullptr) {
                        cout << "Error: No triangle found in locate" << endl;
                        return nullptr;
                    }
                    break;
                }
            }
            if (edge == nullptr) {
                cout << "Error: No edge found in locate" << endl;
                return nullptr;
            }
            current_triangle = next_triangle;
        }
    }
    cout << "Error: No triangle found in locate" << endl;
    return nullptr;
}

template <typename T> void DelaunayTriangulation<T>::triangulate(int n) {
    if (n == -1) {
        n = mesh_vertices.size();
    }
    for (n_cache; n_cache < n; n_cache++) {  
        //Get new vertex
        auto NewVertex = mesh_vertices[n_cache];
        //Locate
        int insert_type;
        auto t = locate(NewVertex, &insert_type);
        //Add interior edges (edges pointing their selves)
        //Add old neighboring triangles to the new edges added
        //Add new edges to the old neighboring triangles
        //Destoy old triangle (with their edges)
        //Legalize edges
        //cout << "Insert type: " << insert_type << endl;
        //Add new vertex to triangulation (add the new 3 or 4 triangles with their empty edges)
        if(insert_type == -1){ // 3 triangles
            auto Pi = t->vertices[0];
            auto Pj = t->vertices[1];
            auto Pk = t->vertices[2];

            auto PiPj = t->edges[0];
            auto PjPk = t->edges[1];
            auto PkPi = t->edges[2];

            //Create new triangles
            auto PiPjNew = new Triangle<T>(Pi, Pj, NewVertex);
            auto PjPkNew = new Triangle<T>(Pj, Pk, NewVertex);
            auto PkPiNew = new Triangle<T>(Pk, Pi, NewVertex);

            //Create new edges
            auto PiNew = new Edge<T>();
            auto PjNew = new Edge<T>();
            auto PkNew = new Edge<T>();

            //Set new edges

            //PiPjNew Triangle
            PiPj->replace_adjacent_triangle(t, TrianglePair(PiPjNew, 0));
            PjNew->add_adjacent_triangle(TrianglePair(PiPjNew, 1));
            PiNew->add_adjacent_triangle(TrianglePair(PiPjNew, 2));

            //PjPkNew Triangle
            PjPk->replace_adjacent_triangle(t, TrianglePair(PjPkNew, 0));
            PkNew->add_adjacent_triangle(TrianglePair(PjPkNew, 1));
            PjNew->add_adjacent_triangle(TrianglePair(PjPkNew, 2));
            PkNew->set_vertices(Pk, NewVertex);

            //PkPiNew Triangle
            PkPi->replace_adjacent_triangle(t, TrianglePair(PkPiNew, 0));
            PiNew->add_adjacent_triangle(TrianglePair(PkPiNew, 1));
            PkNew->add_adjacent_triangle(TrianglePair(PkPiNew, 2));

            //Set vertices
            PiNew->set_vertices(Pi, NewVertex);
            PjNew->set_vertices(Pj, NewVertex);
            PkNew->set_vertices(Pk, NewVertex);

            //Set triangles edges
            PiPjNew->set_edges(PiPj, PjNew, PiNew);
            PjPkNew->set_edges(PjPk, PkNew, PjNew);
            PkPiNew->set_edges(PkPi, PiNew, PkNew);
            PiPjNew->debug_number = n_cache*3;
            PjPkNew->debug_number = n_cache*3 + 1;
            PkPiNew->debug_number = n_cache*3 + 2;

            //add new triangles to mesh
            mesh_triangles.push_back(PiPjNew);
            mesh_triangles.push_back(PjPkNew);
            mesh_triangles.push_back(PkPiNew);

            //add new edges to mesh
            mesh_edges.push_back(PiNew);
            mesh_edges.push_back(PjNew);
            mesh_edges.push_back(PkNew);

            //Destroy old triangle
            mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), t), mesh_triangles.end());
            //free(t);
            t->~Triangle();
            //Legalize
            legalizeEdge(PiPjNew, PiPj, NewVertex);
            legalizeEdge(PjPkNew, PjPk, NewVertex);
            legalizeEdge(PkPiNew, PkPi, NewVertex);       
                


        } else { // 4 triangles or 2 triangles
            auto PjPk = t->edges[insert_type]; //COLINEAR EDGE
            int i_index = (insert_type - 1 + 3) % 3;
            auto Pi = t->vertices[i_index];
            auto Pj = t->vertices[insert_type];
            auto Pk = t->vertices[(insert_type + 1) % 3];
            Triangle<T>* neighbor_triangle = PjPk->getOtherTriangle(t);
            if (neighbor_triangle != nullptr) { // 4 triangles
                Vertex<T>* Pl; // neighbor_triangle->vertices[i]
                int Pl_index;
                for (int i = 0; i < 3; i++) {
                    if(neighbor_triangle->vertices[i] != Pi && neighbor_triangle->vertices[i] != Pj && neighbor_triangle->vertices[i] != Pk){
                        Pl = neighbor_triangle->vertices[i];
                        Pl_index = i;
                        break;
                    }
                }
                auto PiPj = t->edges[(insert_type + 2) % 3];
                //auto PjPk = t->edges[insert_type];
                auto PkPi = t->edges[(insert_type + 1) % 3];
                auto PjPl = neighbor_triangle->edges[(Pl_index - 1 + 3) % 3];
                auto PlPk = neighbor_triangle->edges[Pl_index];

                //create new triangles
                auto PiPjNew = new Triangle<T>(Pi, Pj, NewVertex);
                auto PjPlNew = new Triangle<T>(Pj, Pl, NewVertex);
                auto PlPkNew = new Triangle<T>(Pl, Pk, NewVertex);
                auto PkPiNew = new Triangle<T>(Pk, Pi, NewVertex);

                //create new edges
                auto PiNew = new Edge<T>();
                auto PjNew = new Edge<T>();
                auto PlNew = new Edge<T>();
                auto PkNew = new Edge<T>();

                //set new edges
                //PiPjNew Triangle
                PiPj->replace_adjacent_triangle(t, TrianglePair(PiPjNew, 0));
                PjNew->add_adjacent_triangle(TrianglePair(PiPjNew, 1));
                PiNew->add_adjacent_triangle(TrianglePair(PiPjNew, 2));

                //PjPlNew Triangle
                PjPl->replace_adjacent_triangle(neighbor_triangle, TrianglePair(PjPlNew, 0));
                PlNew->add_adjacent_triangle(TrianglePair(PjPlNew, 1));
                PjNew->add_adjacent_triangle(TrianglePair(PjPlNew, 2));

                //PlPkNew Triangle
                PlPk->replace_adjacent_triangle(neighbor_triangle, TrianglePair(PlPkNew, 0));
                PkNew->add_adjacent_triangle(TrianglePair(PlPkNew, 1));
                PlNew->add_adjacent_triangle(TrianglePair(PlPkNew, 2));

                //PkPiNew Triangle
                PkPi->replace_adjacent_triangle(t, TrianglePair(PkPiNew, 0));
                PiNew->add_adjacent_triangle(TrianglePair(PkPiNew, 1));
                PkNew->add_adjacent_triangle(TrianglePair(PkPiNew, 2));

                //set vertices
                PiNew->set_vertices(Pi, NewVertex);
                PjNew->set_vertices(Pj, NewVertex);
                PlNew->set_vertices(Pl, NewVertex);
                PkNew->set_vertices(Pk, NewVertex);

                //set triangles edges
                PiPjNew->set_edges(PiPj, PjNew, PiNew);
                PjPlNew->set_edges(PjPl, PlNew, PjNew);
                PlPkNew->set_edges(PlPk, PkNew, PlNew);
                PkPiNew->set_edges(PkPi, PiNew, PkNew);
                PiPjNew->debug_number = n_cache*4;
                PjPlNew->debug_number = n_cache*4 + 1;
                PlPkNew->debug_number = n_cache*4 + 2;
                PkPiNew->debug_number = n_cache*4 + 3;

                //add new triangles to mesh
                mesh_triangles.push_back(PiPjNew);
                mesh_triangles.push_back(PjPlNew);
                mesh_triangles.push_back(PlPkNew);
                mesh_triangles.push_back(PkPiNew);

                //add new edges to mesh
                mesh_edges.push_back(PiNew);
                mesh_edges.push_back(PjNew);
                mesh_edges.push_back(PlNew);
                mesh_edges.push_back(PkNew);

                //Destroy old triangles
                mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), t), mesh_triangles.end());
                mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), neighbor_triangle), mesh_triangles.end());
                //free(t);
                //free(neighbor_triangle);
                //t->~Triangle();
                //neighbor_triangle->~Triangle();

                //Destroy old edge
                mesh_edges.erase(std::remove(mesh_edges.begin(), mesh_edges.end(), PjPk), mesh_edges.end());
                //free(PjPk);
                //PjPk->~Edge();
                t->~Triangle();
                neighbor_triangle->~Triangle();
                PjPk->~Edge();
                //Legalize
                legalizeEdge(PiPjNew, PiPj, NewVertex);
                legalizeEdge(PjPlNew, PjPl, NewVertex);
                legalizeEdge(PlPkNew, PlPk, NewVertex);
                legalizeEdge(PkPiNew, PkPi, NewVertex);

            } else { // 2 triangles, but this should not happen for the used boundary triangle
                      
            }
        }
        
    }

    //Save the mesh prio removing the super triangle or inserting constrained edges
    mesh_edges_debug = mesh_edges;
    mesh_triangles_debug = mesh_triangles;
    mesh_vertices_debug = mesh_vertices;
    
}

//Erase boundary triangle and its edges
template<typename T> void DelaunayTriangulation<T>::remove_super_triangle(){
     mesh_triangles_debug.erase(std::remove_if(mesh_triangles_debug.begin(), mesh_triangles_debug.end(), [](Triangle<T>* t) {
        return t->vertices[0]->id < 0 || t->vertices[1]->id < 0 || t->vertices[2]->id < 0;
    }), mesh_triangles_debug.end());

    mesh_edges_debug.erase(std::remove_if(mesh_edges_debug.begin(), mesh_edges_debug.end(), [](Edge<T>* e) {
        return e->adjacent_triangles[0].first->vertices[e->adjacent_triangles[0].second]->id < 0 || 
        e->adjacent_triangles[0].first->vertices[(e->adjacent_triangles[0].second + 1) % 3]->id < 0
        || e->adjacent_triangles[1].first->vertices[e->adjacent_triangles[1].second]->id < 0 || 
        e->adjacent_triangles[1].first->vertices[(e->adjacent_triangles[1].second + 1) % 3]->id < 0;
    }), mesh_edges_debug.end());

    mesh_vertices_debug.erase(std::remove_if(mesh_vertices_debug.begin(), mesh_vertices_debug.end(), [](Vertex<T>* v) {
        return v->id < 0;
    }), mesh_vertices_debug.end());
}

template<typename T> void DelaunayTriangulation<T>::insert_constrained_edges(){

    //Constrained edges
    //Copy of lexico ordered vertices
    std::vector<Vertex<T>*> vertices_copy = mesh_vertices;
    std::sort(vertices_copy.begin(), vertices_copy.end(), [](Vertex<T>* v1, Vertex<T>* v2) {
        if(v1->x_coord == v2->x_coord){
            return v1->y_coord < v2->y_coord;
        }
        return v1->x_coord < v2->x_coord;
    });
    
    //For each constrained edge
    for(auto &cons_e : constraid_edges){
        //Find the colinear vertex with cons_e
        Vertex<T>* ceA = cons_e->vertices[0];
        Vertex<T>* ceB = cons_e->vertices[1];
        cout << "ceA: " << ceA->x_coord << " - " << ceA->y_coord << endl;
        cout << "ceB: " << ceB->x_coord << " - " << ceB->y_coord << endl;
        Predicates<T> pred;
        auto colinear_vertices = pred.find_collinear_points(ceA, ceB, &vertices_copy);
        //cout << "Colinear vertices: " << colinear_vertices.size() << endl;
        //For each pi and pi+1 colinear point
        for(size_t i = 0; i < colinear_vertices.size() - 1; i++){
            auto Pi = colinear_vertices[i];
            auto Pj = colinear_vertices[i+1];
            cout << "Pi: " << Pi->x_coord << " - " << Pi->y_coord << endl;
            cout << "Pj: " << Pj->x_coord << " - " << Pj->y_coord << endl;
            std::vector<std::pair<Triangle<T>*, Edge<T>*>> triangles_to_flip;
            //Find the triangle path
            find_in_path_triangles(&triangles_to_flip, Pi, Pj);
            if(triangles_to_flip.size() == 0){
                cout << "No triangles to flip, segment already exists" << endl;
                continue;
            }
            bool constrained_done = false;
            //Flip the triangles
            cout << "Flipping triangles" << endl;
            cout << "Triangles to and Edges to flip: " << endl;
            for (auto &t : triangles_to_flip) {
                if (t.first != nullptr) {
                    cout << *t.first << endl;
                } else {
                    cout << "Triangle is nullptr" << endl;
                }
                
                if (t.second != nullptr) {
                    cout << *t.second << endl;
                } else {
                    cout << "Edge is nullptr" << endl;
                }
            }
            int debug_counter = 0;
            while(!constrained_done && triangles_to_flip.size() > 0){
                Predicates<T> pred;
                for(int j = 0; j < triangles_to_flip.size()-1; j++){
                    //cout << "Flipping triangle: " << *triangles_to_flip[j].first << endl;
                    //cout << "Flipping edge: " << *triangles_to_flip[j].second << endl;
                    auto t_flip = triangles_to_flip[j].first;
                    auto e_flip = triangles_to_flip[j].second;

                    int legalizing_edge_index = e_flip->get_index_from_triangle(t_flip);//Ojo esto puede ayudar para otras cosas!!!!!!!!!!!!!!!!!
                    if(legalizing_edge_index == -1){
                        cout << "Error: Edge not found in triangle" << endl;
                        break;
                    }

                    if(!e_flip->special){
                        bool flipped_bool = da_way_flip(t_flip, e_flip, {&triangles_to_flip, j}, Pi, Pj);
                        //If the edge was flipped, we need to check if the new edge is PiPj constrained and we have to make it special
                        //in case it is not, this if for disallowing the edge to be flipped again an make the algorithm infinite
                        if(flipped_bool){
                            //Check if new edge is PiPj constrained
                            Edge<T>* PrPk = triangles_to_flip[j].second;
                            Vertex<T>* PrPk_v1 = PrPk->adjacent_triangles[0].first->vertices[PrPk->adjacent_triangles[0].second];
                            Vertex<T>* PrPk_v2 = PrPk->adjacent_triangles[1].first->vertices[PrPk->adjacent_triangles[1].second];
                            if((PrPk_v1 == Pi && PrPk_v2 == Pj) || (PrPk_v1 == Pj && PrPk_v2 == Pi)){
                                constrained_done = true;
                                PrPk->constrained = true;
                            }
                            else PrPk->special = true;

                        }
                        
                    }
                    else{
                        Vertex<T>* P1 = t_flip->vertices[legalizing_edge_index];
                        Vertex<T>* P2 = t_flip->vertices[(legalizing_edge_index + 1) % 3];
                        cout << "Special edge P1: " << *P1 << endl;
                        cout << "Special edge P2: " << *P2 << endl;
                        if(pred.doIntersect_loose(Pi, Pj, P1, P2)){
                            cout << "Special segment was intersecting with PiPj" << endl;
                            e_flip->special = false;
                        }
                    }
                    
                }
                debug_counter++;
                if(debug_counter > 100){
                    cout << "Debug counter break" << endl;
                    break;
                }
            }
        }
    }
    mesh_edges_debug = mesh_edges;
    mesh_triangles_debug = mesh_triangles;
    mesh_vertices_debug = mesh_vertices;
}


template <typename T> void DelaunayTriangulation<T>::find_in_path_triangles(std::vector<std::pair<Triangle<T>*, Edge<T>*>>* triangles_to_flip,
 Vertex<T>* Pi, Vertex<T>* Pj){
    //Find the triangles that contain Pi
    //Linear search, can be optimized
    //We find the first triangle that contains Pi and its neighbors around Pi
    std::vector<Triangle<T>*> triangles;
    //Edge<T>* edge_2;
    //Edge<T>* edge_3;
    for(auto &t : mesh_triangles){
        if(t->contains_vertex(Pi)){
            triangles.push_back(t);
            //break;
        }
    }
    
    //for(int i = 0; i < 3; i++){
    //    if(*(triangles[0]->vertices[i]) == *Pi){
    //        edge_2 = triangles[0]->edges[i];
    //        edge_3 = triangles[0]->edges[(i-1+3)%3];
    //        break;
    //    }
    //}
    //triangles.push_back(edge_2->getOtherTriangle(triangles[0]));
    //triangles.push_back(edge_3->getOtherTriangle(triangles[0]));
    //print triangles
    //From the triangles that contain Pi, find the one that intersects the path PiPj
    for(auto &t : triangles){
        int Pi_index;
        for(int i = 0; i < 3; i++){
            if(*(t->vertices[i]) == *Pi){
                Pi_index = i;
                break;
            }
        }
        Vertex<T>* Pk = t->vertices[(Pi_index + 1) % 3];
        Vertex<T>* Pl = t->vertices[(Pi_index - 1 + 3) % 3];
        Predicates<T> pred;
        Vertex<T> intersection;

        if(t->contains_vertex(Pj)){
            //The path is already created by the triangulation, we extract the edge of the triangle that
            //has PiPj, for that we need to know the orientation of the triangle, or simply search the 
            //Pj index
            int Pj_index;
            for(int i = 0; i < 3; i++){
                if(*(t->vertices[i]) == *Pj){
                    Pj_index = i;
                    break;
                }
            }
            Edge<T>* already_PiPj;
            //PiPj case
            if(Pj_index == (Pi_index + 1) % 3){
                already_PiPj = t->edges[Pi_index];
                
            }
            //PjPi case
            if(Pj_index == (Pi_index - 1 + 3) % 3){
                already_PiPj = t->edges[(Pi_index - 1 + 3) % 3];
                
            }
            already_PiPj->constrained = true;
            cout << "Segment already created by the triangulation, no flips needed" << endl;
            cout << "Edge: " << *already_PiPj << endl;
            cout << "Triangle: " << *t << endl;
            break;
        }
        
        if(pred.doIntersect(Pi, Pj, Pk, Pl, &intersection)){
            triangles_to_flip->push_back(std::make_pair(t, t->edges[(Pi_index+1) % 3]));
            break;
        }

        
        
    }
    cout << "Triangles to flip: " << triangles_to_flip->size() << endl;    
    if(triangles_to_flip->size() == 0){
        cout << "Segment already created by the triangulation, no flips needed" << endl;
        //The edge is already created by the triangulation, we only need to update it with the constrained flag

        return;
    }
    //Auxilary fuction to find the rest of triangle on the track
    cout << "First triangle to flip: " << *triangles_to_flip->back().first << endl;
    this_is_the_way(triangles_to_flip, Pi, Pj);

}

template <typename T> void DelaunayTriangulation<T>::this_is_the_way(std::vector<std::pair<Triangle<T>*, Edge<T>*>>* triangles_to_flip, Vertex<T>* Pi, Vertex<T>* Pj){
    Edge<T>* last_edge = triangles_to_flip->back().second;
    cout << "Last edge: " << *last_edge << endl;
    Triangle<T>* last_triangle = triangles_to_flip->back().first;
    cout << "Last triangle: " << *last_triangle << endl;
    Triangle<T>* current_triangle = last_edge->getOtherTriangle(last_triangle);
    cout << "Current triangle: " << *current_triangle << endl;
    //If Pj is in the triangle, we are done (last triangle to be added, no edge)
    if(current_triangle->contains_vertex(Pj)){
        triangles_to_flip->push_back(std::make_pair(current_triangle, nullptr));
        cout << "End of the path" << endl;
        return;
    }
    //2 edges that can be the next edge, we search for the intersection
    int Pl_index = last_edge->get_index_from_triangle(current_triangle);
    //Notice, we do not care about the original orientation
    std::pair<Vertex<T>*, Vertex<T>*> P2Pl = std::make_pair(current_triangle->vertices[(Pl_index - 1 + 3) % 3], current_triangle->vertices[Pl_index]);
    std::pair<Vertex<T>*, Vertex<T>*> P1P2 = std::make_pair(current_triangle->vertices[(Pl_index + 1) % 3], current_triangle->vertices[(Pl_index - 1 + 3) % 3]);
    //We search for the intersection, the succes case is the next edge recursively
    Predicates<T> pred;
    if(pred.doIntersect(Pi, Pj, P2Pl.first, P2Pl.second)){
        triangles_to_flip->push_back(std::make_pair(current_triangle, current_triangle->edges[(Pl_index - 1 + 3) % 3]));
        this_is_the_way(triangles_to_flip, Pi, Pj);
    } else if(pred.doIntersect(Pi, Pj, P1P2.first, P1P2.second)){
        triangles_to_flip->push_back(std::make_pair(current_triangle, current_triangle->edges[(Pl_index + 1) % 3]));
        this_is_the_way(triangles_to_flip, Pi, Pj);
    } else {
        cout << "Error: No triangle found in this_is_the_way" << endl;
    }

    

}

template <typename T> void DelaunayTriangulation<T>::legalizeEdge(Triangle<T>* PiPjPr, Edge<T>* PiPj, Vertex<T>* v){
    Triangle<T>* PkPjPi = PiPj->getOtherTriangle(PiPjPr);//neighbor triangle
    if (PkPjPi == nullptr) {
        //cout << "Error: Legalizing a boundary edge, can not be legalized" << endl;
        return;
    }
    //Determine the index of the legalizing edge in the original triangle
    int legalizing_edge_index;
    for (int i = 0; i < 3; i++) {
        if(PiPjPr->edges[i] == PiPj){
            legalizing_edge_index = i;
            break;
        }
    }

    //Original triangle vertex
    auto Pi = PiPjPr->vertices[legalizing_edge_index];
    auto Pj = PiPjPr->vertices[(legalizing_edge_index + 1) % 3];
    auto Pr = PiPjPr->vertices[(legalizing_edge_index -1 + 3) % 3];

    //Origninal triangle edges
    //auto PiPj = PiPjPr->edges[0] is the edge we are legalizing
    auto PjPr = PiPjPr->edges[(legalizing_edge_index + 1) % 3];
    auto PrPi = PiPjPr->edges[(legalizing_edge_index -1 + 3) % 3];

    //Neighbor triangle vertex
    Vertex<T>* Pk; // neighbor_triangle->vertices[i]
    int Pk_index;
    for (int i = 0; i < 3; i++) {
        if(PkPjPi->vertices[i] != Pi && PkPjPi->vertices[i] != Pj && PkPjPi->vertices[i] != Pr){
            Pk = PkPjPi->vertices[i];
            Pk_index = i;
            break;
        }
    }
    //Neighbor triangle edges
    auto PiPk = PkPjPi->edges[(Pk_index - 1 + 3) % 3];
    auto PkPj = PkPjPi->edges[Pk_index];

    //Check if the edge is illegal
    Predicates<T> pred;
    if(pred.in_circle(PiPjPr, Pk)){
        //Create new triangles
        auto PiPkPr = new Triangle<T>(Pi, Pk, Pr);
        auto PkPjPr = new Triangle<T>(Pk, Pj, Pr);

        //Create new edge
        auto PrPk = new Edge<T>();

        //Set new edge
        PrPk->add_adjacent_triangle(TrianglePair(PiPkPr, 1));
        PrPk->add_adjacent_triangle(TrianglePair(PkPjPr, 2));

        //Set triangles edges
        PiPkPr->set_edges(PiPk, PrPk, PrPi);
        PkPjPr->set_edges(PkPj, PjPr, PrPk);

        //Set new edge vertices
        PrPk->set_vertices(Pr, Pk);

        //Set old neighboring edges
        PiPk->replace_adjacent_triangle(PkPjPi, TrianglePair(PiPkPr, 0));
        PkPj->replace_adjacent_triangle(PkPjPi, TrianglePair(PkPjPr, 0));
        PjPr->replace_adjacent_triangle(PiPjPr, TrianglePair(PkPjPr, 1));
        PrPi->replace_adjacent_triangle(PiPjPr, TrianglePair(PiPkPr, 2));

        //Remove flipped edge
        mesh_edges.erase(std::remove(mesh_edges.begin(), mesh_edges.end(), PiPj), mesh_edges.end());

        //Remove flipped triangles
        mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), PiPjPr), mesh_triangles.end());
        mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), PkPjPi), mesh_triangles.end());
        //Add new triangles
        mesh_triangles.push_back(PiPkPr);
        mesh_triangles.push_back(PkPjPr);

        //Add new edge
        mesh_edges.push_back(PrPk);

        //Legalize recursively
        legalizeEdge(PiPkPr, PiPk, v);
        legalizeEdge(PkPjPr, PkPj, v);

    }
}


template <typename T>
bool DelaunayTriangulation<T>::da_way_flip(Triangle<T>* PiPjPr, Edge<T>* PiPj, 
    std::pair<std::vector<std::pair<Triangle<T>*, Edge<T>*>>*, int> triangles_to_update, Vertex<T>* Pinit, Vertex<T>* Pfinal) {
        //cout << "Da wae fliping" << endl;
        //cout << "PiPjPr: " << *PiPjPr << endl;
    Predicates<T> pred;
    Triangle<T>* PkPjPi = PiPj->getOtherTriangle(PiPjPr);//neighbor triangle
    if (PkPjPi == nullptr) {
        cout << "Error: Legalizing a boundary edge, can not be legalized" << endl;
        return false;
    }
    
    //cout << "PkPjPi: " << *PkPjPi << endl;
    //Determine the index of the legalizing edge in the original triangle
    int legalizing_edge_index = -1; //Ojo segmentation puede haber aqui, en ese caso hay problema de orientacion
    for (int i = 0; i < 3; i++) {
        if(PiPjPr->edges[i] == PiPj){
            legalizing_edge_index = i;
            break;
        }
    }
    if(legalizing_edge_index == -1){
        cout << "Error: Legalizing edge index is nullptr" << endl;
        return false;
    }

    if(!pred.doIntersect(Pinit, Pfinal, PiPjPr->vertices[(legalizing_edge_index)], PiPjPr->vertices[(legalizing_edge_index +1) % 3])){
        cout << "Segment already created by the triangulation, no flips needed" << endl;
        return true;
    }

    //Original triangle vertex
    auto Pi = PiPjPr->vertices[legalizing_edge_index];
    auto Pj = PiPjPr->vertices[(legalizing_edge_index + 1) % 3];
    auto Pr = PiPjPr->vertices[(legalizing_edge_index -1 + 3) % 3];

    //Origninal triangle edges
    //auto PiPj = PiPjPr->edges[0] is the edge we are legalizing
    auto PjPr = PiPjPr->edges[(legalizing_edge_index + 1) % 3];
    auto PrPi = PiPjPr->edges[(legalizing_edge_index -1 + 3) % 3];

    //Neighbor triangle vertex
    Vertex<T>* Pk; // neighbor_triangle->vertices[i]
    int Pk_index;
    for (int i = 0; i < 3; i++) {
        if(PkPjPi->vertices[i] != Pi && PkPjPi->vertices[i] != Pj && PkPjPi->vertices[i] != Pr){
            Pk = PkPjPi->vertices[i];
            Pk_index = i;
            break;
        }
    }
    //Neighbor triangle edges
    auto PiPk = PkPjPi->edges[(Pk_index - 1 + 3) % 3];
    auto PkPj = PkPjPi->edges[Pk_index];

    //Check if the edge is possible to flips
    //Check if diagonals intersect (if they do, the polygon is convex and the edge is possible)
    Vertex<T> intersection;
    if(pred.doIntersect(Pi, Pj, Pr, Pk, &intersection)){
        //Create new triangles
        auto PiPkPr = new Triangle<T>(Pi, Pk, Pr);
        auto PkPjPr = new Triangle<T>(Pk, Pj, Pr);

        //Create new edge
        auto PrPk = new Edge<T>();

        //Set new edge
        PrPk->add_adjacent_triangle(TrianglePair(PiPkPr, 1));
        PrPk->add_adjacent_triangle(TrianglePair(PkPjPr, 2));

        //Set triangles edges
        PiPkPr->set_edges(PiPk, PrPk, PrPi);
        PkPjPr->set_edges(PkPj, PjPr, PrPk);

        //Set new edge vertices
        PrPk->set_vertices(Pr, Pk);

        //Set old neighboring edges
        PiPk->replace_adjacent_triangle(PkPjPi, TrianglePair(PiPkPr, 0));
        PkPj->replace_adjacent_triangle(PkPjPi, TrianglePair(PkPjPr, 0));
        PjPr->replace_adjacent_triangle(PiPjPr, TrianglePair(PkPjPr, 1));
        PrPi->replace_adjacent_triangle(PiPjPr, TrianglePair(PiPkPr, 2));

        //Remove flipped edge
        mesh_edges.erase(std::remove(mesh_edges.begin(), mesh_edges.end(), PiPj), mesh_edges.end());

        //Remove flipped triangles
        mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), PiPjPr), mesh_triangles.end());
        mesh_triangles.erase(std::remove(mesh_triangles.begin(), mesh_triangles.end(), PkPjPi), mesh_triangles.end());

        //Add new triangles
        mesh_triangles.push_back(PiPkPr);
        mesh_triangles.push_back(PkPjPr);

        //Add new edge
        mesh_edges.push_back(PrPk);
        Triangle<T>* next_triangle = PiPkPr;
        Triangle<T>* past_triangle = PkPjPr;
        Edge<T>* next_flip_edge = nullptr;
        
        if((*triangles_to_update.first)[triangles_to_update.second + 1].second != nullptr){
            //Identify the first triangle and the second
            next_flip_edge = (*triangles_to_update.first)[triangles_to_update.second + 1].second;
            for(int i = 0; i < 3; i++){
                if(past_triangle->edges[i] == next_flip_edge){
                    next_triangle = PkPjPr;
                    past_triangle = PiPkPr;
                    
                }
            }
        }

        //Triangles flipped, so the first one is not in the way anymore
        //past_triangle->in_da_wae = false;
        //cout << *past_triangle << endl;
        //cout << *next_triangle << endl;

        //Update the triangles to flip vector
        (*triangles_to_update.first)[triangles_to_update.second] = std::make_pair(past_triangle, PrPk);
        
        //The next edge keeps the same, we just changed the i and i+1 triangles and the i edge
        (*triangles_to_update.first)[triangles_to_update.second + 1] = std::make_pair(next_triangle, next_flip_edge);

        //Print all triangles
        //cout << "Triangles and Edges after flip" << endl;
        //cout << endl;
        //print_triangles();
        //print_edges();

        return true;

    }
    return false;
    
}


template <typename T> void DelaunayTriangulation<T>::print_input_vertex() {
    cout << "Input Vertex: " << endl;
    for (auto v : mesh_vertices) {
        std::cout << *v << std::endl;
    }
}

template <typename T> void DelaunayTriangulation<T>::print_edges() {
    cout << "Edges: " << endl;
    for (auto e : mesh_edges) {
        std::cout << *e << std::endl;
        if (e->adjacent_triangles[0].first != nullptr && e->adjacent_triangles[1].first != nullptr){
            std::cout << e->adjacent_triangles[0].first->debug_number << " " << e->adjacent_triangles[1].first->debug_number << std::endl;
        }
        else if (e->adjacent_triangles[0].first != nullptr){
            std::cout << e->adjacent_triangles[0].first->debug_number << " " << "NULL" << std::endl;
        }
        else if (e->adjacent_triangles[1].first != nullptr){
            std::cout << "NULL" << " " << e->adjacent_triangles[1].first->debug_number << std::endl;
        }
        else{
            std::cout << "NULL" << " " << "NULL" << std::endl;
        }
    }
}

template <typename T> void DelaunayTriangulation<T>::print_triangles() {
    cout << "Triangles: " << endl;
    for (auto t : mesh_triangles) {
        std::cout << "ID: "<< t->debug_number << *t << std::endl;
    }
}

template <typename T> std::vector<Triangle<T>> DelaunayTriangulation<T>::getTriangles() {
    std::vector<Triangle<T>> triangles;
    for (auto t : mesh_triangles) {
        triangles.push_back(*t);
    }
    return triangles;
}

template <typename T> std::vector<Edge<T>> DelaunayTriangulation<T>::getEdges() {
    std::vector<Edge<T>> edges;
    for (auto e : mesh_edges) {
        edges.push_back(*e);
    }
    return edges;
}

template <typename T> std::vector<Vertex<T>> DelaunayTriangulation<T>::getVertices() {
    std::vector<Vertex<T>> vertices;
    for (auto v : mesh_vertices) {
        vertices.push_back(*v);
    }
    return vertices;
}

#endif // DELAUNAY_TRIANGULATION_HPP