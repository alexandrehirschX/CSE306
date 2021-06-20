#pragma once

#include "vector.hpp"

class Polygon{
    public:
        std::vector<Vector> vertices;
        std::vector<std::pair<Vector,Vector>> edges;

    explicit Polygon(std::vector<Vector> V = {}){
        vertices = V;
        for (size_t i = 0; i < vertices.size(); i++){
            // std::pair<Vector, Vector> edge(V[i],V[(i != vertices.size()-1) ? i+1 : 0]);
            std::pair<Vector, Vector> edge(V[i], V[(i != 0) ? i - 1 : vertices.size() - 1]);
            edges.push_back(edge);
        }
    }
};

Vector intersect(Vector prevVertex, Vector curVertex, std::pair<Vector, Vector> clipEddge);

Vector intersect_voronoi(Vector &prevVertex, Vector &curVertex, Vector &curPoint, Vector &tempPoint, Vector M);

bool inside(Vector P, std::pair<Vector, Vector> clipEdge);

bool inside(Vector X, Vector P, Vector Pj, Vector M);

Polygon clip_poly(Polygon subjectPolygon, Polygon clipPolygon);

std::vector<Polygon> voronoi(std::vector<Vector> points, std::vector<double> weight = {}, Polygon space = Polygon({Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)}));

std::vector<Vector> voronoi_tesselation(std::vector<Vector> points, int Iter);

void save_svg(std::string filename, const std::vector<Polygon> &polygons = {}, const std::vector<Vector> &vectors = {}, std::string fillcol = "none");

void save_svg(std::string filename, const std::vector<Vector> &vectors = {}, const std::vector<Polygon> &polygons = {}, std::string fillcol = "none");

void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes);