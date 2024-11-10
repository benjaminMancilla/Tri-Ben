#include "delaunay-triangulation.hpp"
#include <SFML/Graphics.hpp>

//random points inside a 1x1 square centred in 0.0 using time as seed
std::vector<Vertex<double>> random_points(int n, double radius = 1.0) {
    std::vector<Vertex<double>> points;
    points.reserve(n);  // Reservar espacio para evitar realocaciones

    // Semilla aleatoria basada en el tiempo
    srand(time(NULL));

    for (int i = 0; i < n; ++i) {
        // Generar un ángulo aleatorio entre 0 y 2π
        double angle = ((double)rand() / RAND_MAX) * 2 * M_PI;
        
        // Generar una distancia aleatoria desde el centro hasta el borde del círculo
        // Para que la distribución sea uniforme en el área, usamos la raíz cuadrada
        double distance = radius * sqrt((double)rand() / RAND_MAX);
        
        // Convertir de coordenadas polares a cartesianas
        double x = distance * cos(angle);
        double y = distance * sin(angle);

        // Añadir el punto a la lista
        points.push_back(Vertex<double>(x, y));
    }

    return points;
}

std::vector<Vertex<double>> random_points_in_triangle(double large, int n) {
    std::vector<Vertex<double>> points;
    points.reserve(n + 4);  // Incluir espacio para vértices

    // Definir los cuatro vértices del cuadrilátero centrado en el origen
    points.push_back(Vertex<double>(0, large/sqrt(3)));  // Vértice superior derecho
    points.push_back(Vertex<double>(-large / 2, -large / (2 * sqrt(3))));  // Vértice inferior derecho
    points.push_back(Vertex<double>(large / 2, -large / (2 * sqrt(3))));  // Vértice inferior izquierdo


    // Semilla aleatoria basada en el tiempo
    srand(static_cast<unsigned int>(time(NULL)));
    double radius = large * 2 / sqrt(3);

    for (int i = 0; i < n; ++i) {
        // Generar un ángulo aleatorio entre 0 y 2π
        double angle = ((double)rand() / RAND_MAX) * 2 * M_PI;
        
        // Generar una distancia aleatoria desde el centro hasta el borde del círculo
        // Para que la distribución sea uniforme en el área, usamos la raíz cuadrada
        double distance = radius * sqrt((double)rand() / RAND_MAX);
        
        // Convertir de coordenadas polares a cartesianas
        double x = distance * cos(angle);
        double y = distance * sin(angle);

        // Añadir el punto a la lista
        points.push_back(Vertex<double>(x, y));
    }

    return points;
}

std::vector<Vertex<double>> random_points_in_quadrilateral(double width, double height, int n) {
    std::vector<Vertex<double>> points;
    points.reserve(n + 4);  // Incluir espacio para vértices

    // Definir los cuatro vértices del cuadrilátero centrado en el origen
    points.push_back(Vertex<double>( width / 2,  height / 2));  // Vértice superior derecho
    points.push_back(Vertex<double>( width / 2, -height / 2));  // Vértice inferior derecho
    points.push_back(Vertex<double>(-width / 2, -height / 2));  // Vértice inferior izquierdo
    points.push_back(Vertex<double>(-width / 2,  height / 2));  // Vértice superior izquierdo

    // Semilla aleatoria basada en el tiempo
    srand(static_cast<unsigned int>(time(NULL)));

    // Generar n puntos aleatorios dentro del cuadrilátero
    for (int i = 0; i < n; ++i) {
        // Generar una posición aleatoria dentro de los límites del cuadrilátero
        double x = (static_cast<double>(rand()) / RAND_MAX) * width - width / 2;
        double y = (static_cast<double>(rand()) / RAND_MAX) * height - height / 2;

        points.push_back(Vertex<double>(x, y));
    }

    return points;
}

std::vector<Vertex<double>> random_points_in_C(double large, double thick, int n) {
    std::vector<Vertex<double>> points;
    points.reserve(n + 4);  // Incluir espacio para vértices

    // Definir los cuatro vértices del cuadrilátero centrado en el origen
    points.push_back(Vertex<double>( large / 2,  large / 2));  // Vértice superior derecho
    points.push_back(Vertex<double>(large / 2,  (large / 2 - thick)));  // Vértice superior derecho
    points.push_back(Vertex<double>((-large / 2 + thick),  (large / 2 - thick)));  // Vértice superior izquierdo
    points.push_back(Vertex<double>((-large / 2 + thick),  (-large / 2 + thick)));  // Vértice inferior izquierdo
    points.push_back(Vertex<double>(large / 2,  (-large / 2 + thick)));  // Vértice inferior derecho
    points.push_back(Vertex<double>( large / 2, -large / 2));  // Vértice inferior derecho
    points.push_back(Vertex<double>(-large / 2, -large / 2));  // Vértice inferior izquierdo
    points.push_back(Vertex<double>(-large / 2,  large / 2));  // Vértice superior izquierdo
   
    

    // Semilla aleatoria basada en el tiempo
    srand(static_cast<unsigned int>(time(NULL)));
    double radius = large * 2;

    // Generar n puntos aleatorios dentro del cuadrilátero
    for (int i = 0; i < n; ++i) {
        // Generar un ángulo aleatorio entre 0 y 2π
        double angle = ((double)rand() / RAND_MAX) * 2 * M_PI;
        
        // Generar una distancia aleatoria desde el centro hasta el borde del círculo
        // Para que la distribución sea uniforme en el área, usamos la raíz cuadrada
        double distance = radius * sqrt((double)rand() / RAND_MAX);
        
        // Convertir de coordenadas polares a cartesianas
        double x = distance * cos(angle);
        double y = distance * sin(angle);

        // Añadir el punto a la lista
        points.push_back(Vertex<double>(x, y));
    }

    return points;
}

std::vector<Vertex<double>> random_points_in_quadrilateral_with_line(double width, double height, int n) {
    std::vector<Vertex<double>> points;
    points.reserve(n + 4);  // Incluir espacio para vértices

    // Definir los cuatro vértices del cuadrilátero centrado en el origen
    points.push_back(Vertex<double>( width / 2,  height / 2));  // Vértice superior derecho
    points.push_back(Vertex<double>( width / 2, -height / 2));  // Vértice inferior derecho
    points.push_back(Vertex<double>(-width / 2, -height / 2));  // Vértice inferior izquierdo
    points.push_back(Vertex<double>(-width / 2,  height / 2));  // Vértice superior izquierdo
    points.push_back(Vertex<double>(width / 4, 0));  // Vértice en el centro
    points.push_back(Vertex<double>(-width / 4, 0));  // Vértice en el centro superior

    // Semilla aleatoria basada en el tiempo
    srand(static_cast<unsigned int>(time(NULL)));

    // Generar n puntos aleatorios dentro del cuadrilátero
    for (int i = 0; i < n; ++i) {
        // Generar una posición aleatoria dentro de los límites del cuadrilátero
        double x = (static_cast<double>(rand()) / RAND_MAX) * width - width / 2;
        double y = (static_cast<double>(rand()) / RAND_MAX) * height - height / 2;

        points.push_back(Vertex<double>(x, y));
    }

    return points;
}


std::vector<Vertex<double>> grid_points(int num_points_per_side, double spacing = 1.0) {
    std::vector<Vertex<double>> points;
    points.reserve(num_points_per_side * num_points_per_side);  // Reservar espacio

    double half_size = (num_points_per_side - 1) * spacing / 2.0;  // Mitad del tamaño total de la grilla

    for (int i = 0; i < num_points_per_side; ++i) {
        for (int j = 0; j < num_points_per_side; ++j) {
            // Las coordenadas X y Y del punto, centrando la grilla en el origen (0,0)
            double x = i * spacing - half_size;
            double y = j * spacing - half_size;

            points.push_back(Vertex<double>(x, y));
        }
    }

    return points;
}

std::vector<Vertex<double>> colinear_test_points() {
    std::vector<Vertex<double>> points;
    points.push_back(Vertex<double>(-1, -1));
    points.push_back(Vertex<double>(1, 1));
    points.push_back(Vertex<double>(-1, 1));
    points.push_back(Vertex<double>(1, -1));
    points.push_back(Vertex<double>(0.5, 0.5));
    points.push_back(Vertex<double>(-0.5, -0.5));
    return points;
}

void display_triangulation(DelaunayTriangulation<double> &triangulation) {
    sf::RenderWindow window(sf::VideoMode(800, 800), "Delaunay Triangulation");
    window.setFramerateLimit(60);

    // Cargar una fuente para los números en los ejes
    sf::Font font;
    if (!font.loadFromFile("arial.ttf")) {
        // Manejo de error si no se encuentra la fuente
        return;
    }

    // Configuración inicial
    window.clear(sf::Color::Black);
    sf::View view = window.getView();
    float zoomFactor = 1.0f;

    // Evento de zoom
    sf::Event event;
    bool running = true;

    std::unordered_set<Triangle<double>*> LEPP_set;
    Edge<double>* terminal_edge_render = nullptr;

    while (running) {
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                running = false;
                window.close();
            }

            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                // Convertir coordenadas de pixel a coordenadas del sistema de la vista actual
                sf::Vector2f clickPos = window.mapPixelToCoords(sf::Mouse::getPosition(window), view);
                
                // Aplicar el factor de escala inverso usado en la visualización y ajustar eje Y
                double scaledX = clickPos.x / 80.0;
                double scaledY = -clickPos.y / 80.0;  // Negar el Y para invertir el eje
                Vertex<double> clickVertex(scaledX, scaledY);

                std::cout << "clickVertex: " << clickVertex << std::endl;
                
                std::cout << "Click en coordenadas de triangulación: (" << scaledX << ", " << scaledY << ")" << std::endl;

                // Buscamos el Triángulo que contiene el click
                Triangle<double>* t0 = triangulation.locate(&clickVertex);

                if (t0 != nullptr) {
                    std::cout << "Triángulo T0 LEPP: " << *t0 << std::endl;
                    // Aquí se puede continuar con LEPP si es necesario
                } else {
                    std::cout << "Triángulo no encontrado" << std::endl;
                }

                //Buscamos el LEPP para el triangulo encontrado
                auto [LEPP_path, terminal_edge] = triangulation.LEPP_search(t0);

                // Agregar los triángulos al conjunto para resaltado
                LEPP_set.clear();
                for (auto t : LEPP_path) {
                    LEPP_set.insert(t);
                }         

                terminal_edge_render = terminal_edge;

                std::cout << "LEPP Path: " << std::endl;
                for (auto t : LEPP_path) {
                    std::cout << *t << std::endl;
                }
  
            }

            // Zoom in ('z') y zoom out ('x')
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Z) {
                    zoomFactor *= 0.9f;  // Zoom in
                    view.zoom(0.9f);
                } else if (event.key.code == sf::Keyboard::X) {
                    zoomFactor *= 1.1f;  // Zoom out
                    view.zoom(1.1f);
                } else if (event.key.code == sf::Keyboard::Escape) {
                    running = false;
                    window.close();
                } else if (event.key.code == sf::Keyboard::Right) {
                    view.move(10, 0);
                } else if (event.key.code == sf::Keyboard::Left) {
                    view.move(-10, 0);
                } else if (event.key.code == sf::Keyboard::Up) {
                    view.move(0, -10);
                } else if (event.key.code == sf::Keyboard::Down) {
                    view.move(0, 10);
                }
            }
        }

        // Redibujar la ventana
        window.clear(sf::Color::Black);
        window.setView(view);

        // Dibujar la grilla (con líneas cada 5 unidades)
        for (int i = -100; i <= 100; i += 5) {
            // Eje vertical
            sf::Vertex verticalLine[] = {
                sf::Vertex(sf::Vector2f(i * 80, -8000), sf::Color(150, 150, 150)),
                sf::Vertex(sf::Vector2f(i * 80, 8000), sf::Color(150, 150, 150))
            };
            //window.draw(verticalLine, 2, sf::Lines);

            // Eje horizontal
            sf::Vertex horizontalLine[] = {
                sf::Vertex(sf::Vector2f(-8000, i * 80), sf::Color(150, 150, 150)),
                sf::Vertex(sf::Vector2f(8000, i * 80), sf::Color(150, 150, 150))
            };
            //window.draw(horizontalLine, 2, sf::Lines);

            // Mostrar los números en los ejes
            if (i % 5 == 0) {
                // Eje x
                sf::Text textX;
                textX.setFont(font);
                textX.setString(std::to_string(i));
                textX.setCharacterSize(15);
                textX.setFillColor(sf::Color::White);
                textX.setPosition(i * 80, 0);  // Ajusta la posición del texto
                window.draw(textX);

                // Eje y
                sf::Text textY;
                textY.setFont(font);
                textY.setString(std::to_string(i));
                textY.setCharacterSize(15);
                textY.setFillColor(sf::Color::White);
                textY.setPosition(0, i * 80);  // Ajusta la posición del texto
                window.draw(textY);
            }
        }

        // Marcar el origen (centro de coordenadas)
        sf::CircleShape originMarker(5);
        originMarker.setFillColor(sf::Color::Yellow);
        originMarker.setPosition(-5, -5);  // Centramos el marcador
        window.draw(originMarker);

        // Dibujar los puntos
        for (const auto &vertex : triangulation.mesh_vertices_debug) {
            sf::CircleShape point(5);
            point.setFillColor(sf::Color::Red);
            point.setPosition(vertex->x_coord * 80, -vertex->y_coord * 80); // Invertir el eje Y
            point.setOrigin(5, 5);  // Centramos el punto
            window.draw(point);
        }

        for (const auto &triangle : triangulation.mesh_triangles_debug) {

            sf::Color triangleColor = sf::Color::White;

            // Dibujar el triángulo
            for (int i = 0; i < 3; i++) {
                Vertex<double>* first_vertex = triangle->vertices[i];
                Vertex<double>* second_vertex = triangle->vertices[(i + 1) % 3];

                sf::Vertex line[] = {
                    sf::Vertex(sf::Vector2f(first_vertex->x_coord * 80, -first_vertex->y_coord * 80), triangleColor),
                    sf::Vertex(sf::Vector2f(second_vertex->x_coord * 80, -second_vertex->y_coord * 80), triangleColor)
                };
                window.draw(line, 2, sf::Lines);
            }
        }

        for (const auto &triangle : LEPP_set) {
            sf::Color triangleColor = sf::Color::Cyan;

            // Dibujar el triángulo
            for (int i = 0; i < 3; i++) {
                Vertex<double>* first_vertex = triangle->vertices[i];
                Vertex<double>* second_vertex = triangle->vertices[(i + 1) % 3];

                sf::Vertex line[] = {
                    sf::Vertex(sf::Vector2f(first_vertex->x_coord * 80, -first_vertex->y_coord * 80), triangleColor),
                    sf::Vertex(sf::Vector2f(second_vertex->x_coord * 80, -second_vertex->y_coord * 80), triangleColor)
                };
                window.draw(line, 2, sf::Lines);
            }
        }

        for (const auto &edge : triangulation.mesh_edges_debug) {
            bool isTerminalEdge = (edge == terminal_edge_render);
            if(edge->constrained || isTerminalEdge){
                
                Vertex<double>* first_vertex = edge->vertices[0];
                Vertex<double>* second_vertex = edge->vertices[1];

                sf::Color edgeColor = isTerminalEdge ? sf::Color::Red : (edge->constrained ? sf::Color::Green : sf::Color::White);

                sf::Vertex line[] = {
                    sf::Vertex(sf::Vector2f(first_vertex->x_coord * 80, -first_vertex->y_coord * 80), edgeColor),
                    sf::Vertex(sf::Vector2f(second_vertex->x_coord * 80, -second_vertex->y_coord * 80), edgeColor)
                };
                window.draw(line, 2, sf::Lines);


            }
        }


        // Mostrar todo lo dibujado
        window.display();
    }
}


void display_points(std::vector<Vertex<double>> &points) {
    sf::RenderWindow window(sf::VideoMode(800, 800), "Delaunay Triangulation");
    window.setFramerateLimit(60);

    // Cargar una fuente para los números en los ejes
    sf::Font font;
    if (!font.loadFromFile("arial.ttf")) {
        // Manejo de error si no se encuentra la fuente
        return;
    }

    // Configuración inicial
    window.clear(sf::Color::Black);
    sf::View view = window.getView();
    float zoomFactor = 1.0f;

    // Evento de zoom
    sf::Event event;
    bool running = true;

    while (running) {
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                running = false;
                window.close();
            }
            // Zoom in ('z') y zoom out ('x')
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Z) {
                    zoomFactor *= 0.9f;  // Zoom in
                    view.zoom(0.9f);
                } else if (event.key.code == sf::Keyboard::X) {
                    zoomFactor *= 1.1f;  // Zoom out
                    view.zoom(1.1f);
                } else if (event.key.code == sf::Keyboard::Escape) {
                    running = false;
                    window.close();
                } else if (event.key.code == sf::Keyboard::Right)
                {
                    view.move(10, 0);
                }
                else if (event.key.code == sf::Keyboard::Left)
                {
                    view.move(-10, 0);
                }
                else if (event.key.code == sf::Keyboard::Up)
                {
                    view.move(0, -10);
                }
                else if (event.key.code == sf::Keyboard::Down)
                {
                    view.move(0, 10);
                }
            }
        }

        // Redibujar la ventana
        window.clear(sf::Color::Black);
        window.setView(view);

        // Dibujar la grilla (con líneas cada 5 unidades)
        for (int i = -100; i <= 100; i += 5) {
            // Eje vertical
            sf::Vertex verticalLine[] = {
                sf::Vertex(sf::Vector2f(i * 80, -8000), sf::Color(150, 150, 150, 255)),
                sf::Vertex(sf::Vector2f(i * 80, 8000), sf::Color(150, 150, 150, 255))
            };
            //window.draw(verticalLine, 2, sf::Lines);

            // Eje horizontal
            sf::Vertex horizontalLine[] = {
                sf::Vertex(sf::Vector2f(-8000, i * 80), sf::Color(150, 150, 150)),
                sf::Vertex(sf::Vector2f(8000, i * 80), sf::Color(150, 150, 150))
            };
            //window.draw(horizontalLine, 2, sf::Lines);

            // Mostrar los números en los ejes
            /*
            if (i % 5 == 0) {
                // Eje x
                sf::Text textX;
                textX.setFont(font);
                textX.setString(std::to_string(i));
                textX.setCharacterSize(15);
                textX.setFillColor(sf::Color::White);
                textX.setPosition(i * 80, 0);  // Ajusta la posición del texto
                window.draw(textX);

                // Eje y
                sf::Text textY;
                textY.setFont(font);
                textY.setString(std::to_string(i));
                textY.setCharacterSize(15);
                textY.setFillColor(sf::Color::White);
                textY.setPosition(0, i * 80);  // Ajusta la posición del texto
                window.draw(textY);
            }
            */
        }

        // Marcar el origen (centro de coordenadas)
        sf::CircleShape originMarker(5);
        originMarker.setFillColor(sf::Color::Yellow);
        originMarker.setPosition(-5, -5);  // Centramos el marcador
        window.draw(originMarker);

        // Dibujar los puntos
        for (auto &p : points) {
            sf::CircleShape point(0.04);
            point.setFillColor(sf::Color::Red);
            point.setPosition(p.x_coord * 80, -p.y_coord * 80); // Invertir el eje Y
            point.setOrigin(5, 5);  // Centramos el punto
            window.draw(point);
        }

        // Mostrar todo lo dibujado
        window.display();
    }
}

void consecutive_constrained_edges(DelaunayTriangulation<double> &dt, int n, int m) {
    std::vector<Vertex<double>> points = random_points(n, 100);
    std::random_shuffle(points.begin(), points.end());
    std::vector<Edge<double>> edges;
    dt.initVertex_and_ConsEdges(points, edges);
    dt.boundingTriangle();
    dt.triangulate();
    for (int i = 0; i < m; i++) {
        dt.constraid_edges.clear();
        Vertex<double> *v1 = dt.mesh_vertices_debug[i];
        Vertex<double> *v2 = dt.mesh_vertices_debug[i + m + 1];
        dt.add_constrainedEdge(v1, v2);
        dt.insert_constrained_edges();
        dt.remove_super_triangle();
        display_triangulation(dt);
        

        
    }
}


int main() {
    DelaunayTriangulation<double> dt;
    int number_of_points = 100;
    int num_points_per_side = 6;  // Grilla de 10x10 puntos
    double spacing = 1.0;  // Distancia entre puntos vecinos

    //std::vector<Vertex<double>> grid = grid_points(num_points_per_side, spacing);
    //std::vector<Vertex<double>> points = random_points(number_of_points, 10);
    //std::vector<Vertex<double>> points = random_points_in_quadrilateral(20, 10, number_of_points/2);
    std::vector<Vertex<double>> points = random_points_in_quadrilateral_with_line(20, 10, number_of_points/2);
    //std::vector<Vertex<double>> points2 = random_points(20, 30);
    //std::vector<Vertex<double>> points = random_points_in_C(20, 3, number_of_points);
    //std::vector<Vertex<double>> points = random_points_in_triangle(20, number_of_points);
    Edge<double> e1 = Edge<double>(&points[0], &points[1]);
    Edge<double> e2 = Edge<double>(&points[1], &points[2]);
    Edge<double> e3 = Edge<double>(&points[2], &points[3]);
    Edge<double> e4 = Edge<double>(&points[3], &points[0]);

    

    std::vector<Edge<double>> edges;
    edges.push_back(e1);
    edges.push_back(e2);
    edges.push_back(e3);
    edges.push_back(e4);

    std::vector<Vertex<double>> all_points = points;
    //all_points.insert(all_points.end(), points2.begin(), points2.end());
    //random order for all points
    std::random_shuffle(all_points.begin(), all_points.end());
    //edges.push_back(e1);    
    dt.initVertex_and_ConsEdges(all_points, edges);

    dt.boundingTriangle();
    //dt.print_input_vertex();
    
    dt.triangulate();
    //display_triangulation(dt);
    dt.insert_constrained_edges();
    dt.remove_super_triangle();
    display_triangulation(dt);
    
    //display_triangulation(dt);
    //dt.remove_super_triangle();
    //cout << "Triangulation done" << endl;
    //dt.print_edges();
    //dt.print_triangles();

   /*
    P1: (-10, -5)
    Q1: (-10, 5)
    P2: (-13.7073, 1.6922)
    Q2: (-9.52029, -1.95206)

    */

    /*
    Vertex<double> P1(-10, -5);
    Vertex<double> Q1(-10, 5);
    Vertex<double> P2(-13.7073, 1.6922);
    Vertex<double> Q2(-9.52029, -1.95206);
    Predicates<double> pred;
    Vertex<double> intersection;
    bool inter = pred.doIntersect(&P1, &Q1, &P2, &Q2, &intersection);
    cout << "Intersection: " << inter << endl;
    cout << "Intersection point: " << intersection << endl;
    */
    
    
    //consecutive_constrained_edges(dt, 100, 3);


    
    //display_triangulation(dt);







    




    return 0;
}