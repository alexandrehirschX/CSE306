#include "object.hpp"

#include <map>
#include <algorithm>

Vector intersect(Vector prevVertex, Vector curVertex, std::pair<Vector,Vector> clipEdge){
    Vector u = clipEdge.first;
    Vector v = clipEdge.second;
    Vector N = Vector(v[1]-u[1],u[0]-v[0],0);
    double t = dot(u-prevVertex,N)/dot(curVertex-prevVertex,N);
    if (t < 0 || t > 1){
        return Vector(inf,inf,inf);
    }
    return prevVertex + t*(curVertex-prevVertex);
};

Vector intersect_voronoi(Vector &A, Vector &B, Vector &Pi, Vector &Pj, Vector M){
    // Vector M = (Pi + Pj)/2;
    double t = dot(M-A,Pi-Pj)/dot(B-A,Pi-Pj);
    return A + t*(B-A);
}

bool inside(Vector P, std::pair<Vector,Vector> clipEdge){
    Vector u = clipEdge.first;
    Vector v = clipEdge.second;
    Vector N = Vector(v[1] - u[1], u[0] - v[0],0);
    return dot(P-u,N) <= 0 ? true : false;
}

bool inside(Vector X, Vector Pi, Vector Pj, Vector M){
    return dot(X-M,Pj-Pi)<0;
}

double max_dist(Polygon polygon,Vector P){
    double dist = 0;
    for (auto &vert: polygon.vertices){
        dist = std::max(dist, distance_square(vert, P));
    }
    return dist;
}

Polygon clip_poly(Polygon subjectPolygon, Polygon clipPolygon){
    Polygon outPolygon;    
    for(auto &clipEdge : clipPolygon.edges){
        outPolygon = Polygon();
        for (size_t i = 0; i < subjectPolygon.vertices.size(); i++){
            Vector curVertex = subjectPolygon.vertices[i];
            Vector prevVertex = subjectPolygon.vertices[(i>0)?(i-1):subjectPolygon.vertices.size()-1];
            Vector intersection = intersect(prevVertex,curVertex,clipEdge);
            if (inside(curVertex,clipEdge)){
                if (!inside(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
                outPolygon.vertices.push_back(curVertex);
            }
            else if (inside(prevVertex,clipEdge)) outPolygon.vertices.push_back(intersection);
        }
        subjectPolygon = outPolygon;
    }
    return outPolygon;
};

std::vector<Polygon> voronoi(std::vector<Vector> points, std::vector<double> weight, Polygon space){
    std::vector<Polygon> res;
    Polygon outPolygon;
    if (weight.size() == 0){
        for (int i = 0; i < points.size(); i++)
        {
            weight.push_back(0.1);
        }
    }
    // #pragma omp parallel for schedule(static,1)
    // for(auto &curPoint: points){
    for (int i = 0; i < points.size();i++){
        auto curPoint = points[i]; 
        std::map<double,Vector> otherPoints;
        std::vector<Vector> otherPoints_vector;
        for(auto &tempPoint: points){
            if (tempPoint != curPoint){
                otherPoints[distance_square(curPoint,tempPoint)] = tempPoint;
                otherPoints_vector.push_back(tempPoint);
            }
        }
        outPolygon = space;
        // int j = 0;
        // for(auto &it: otherPoints){
            // if (it.first > 4*max_dist(outPolygon,curPoint)) break;
                                                                        
            // auto tempPoint = it.second;
        for (int j = 0; j < otherPoints_vector.size();j++){
            Vector tempPoint = otherPoints_vector[j];
            size_t n = outPolygon.vertices.size();
            Polygon tempPolygon = Polygon();
            for(size_t t = 0; t < n; t++){
                Vector curVertex = outPolygon.vertices[t];
                Vector prevVertex = outPolygon.vertices[(t>0)?(t-1):n-1];
                Vector M_temp = (curPoint + tempPoint)/2;
                Vector M = M_temp + (weight[i] - weight[j])/(2*norm_square(curPoint - tempPoint)) * (tempPoint - curPoint);
                Vector intersection = intersect_voronoi(prevVertex,curVertex,curPoint,tempPoint,M);
                if(inside(curVertex,curPoint,tempPoint,M)){
                    if (!inside(prevVertex, curPoint,tempPoint, M))
                        tempPolygon.vertices.push_back(intersection);
                    tempPolygon.vertices.push_back(curVertex);
                }
                else if (inside(prevVertex, curPoint,tempPoint, M))
                    tempPolygon.vertices.push_back(intersection);
            }
            outPolygon = tempPolygon;
            // j += 1;
        }
        res.push_back(outPolygon);
    }
    return res;
}

// std::vector<Polygon> voronoi(std::vector<Vector> points,Polygon space){
//     std::vector<Polygon> res;
//     Polygon outPolygon;
//     for(auto &curPoint: points){
//         outPolygon = space;
//         for(auto &tempPoint: points){
//             if (tempPoint != curPoint){
//                 Vector M = (tempPoint + curPoint)/2;
//                 size_t n = outPolygon.vertices.size();
//                 Polygon tempPolygon = Polygon();
//                 for(int i = 0; i < n; i++){
//                     Vector curVertex = outPolygon.vertices[i];
//                     Vector prevVertex = outPolygon.vertices[(i>0)?(i-1):n-1];
//                     Vector intersection = intersect_voronoi(prevVertex,curVertex,curPoint,tempPoint);
//                     if(inside(curVertex,curPoint,tempPoint)){
//                         if (!inside(prevVertex, curPoint,tempPoint))
//                             tempPolygon.vertices.push_back(intersection);
//                         tempPolygon.vertices.push_back(curVertex);
//                     }
//                     else if (inside(prevVertex, curPoint,tempPoint))
//                         tempPolygon.vertices.push_back(intersection);
//                 }
//                 outPolygon = tempPolygon;
//             }
//         }
//         res.push_back(outPolygon);
//     }
//     return res;
// }


std::vector<Vector> voronoi_tesselation(std::vector<Vector> points, int iter = 20){
    if (iter == 0) return points;
    std::vector<Polygon> vor = voronoi(points);
    std::vector<Vector> new_points;
    for (auto &cell: vor){
        int n = cell.vertices.size();
        double A = 0;   
        Vector C = Vector(); 
        for (size_t i = 0; i < n;i++){
            Vector cur = cell.vertices[i];
            auto next = cell.vertices[(i == n-1) ? 0 : (i+1)];
            A += (cur[0]*next[1] - next[0]*cur[1]);
        }
        A /= 2;
        for (size_t i = 0; i < n;i++){
            auto cur = cell.vertices[i];
            auto next = cell.vertices[(i == n-1) ? 0 : (i+1)];
            C[0] += (cur[0]+next[0])*(cur[0]*next[1]-next[0]*cur[1]);
            C[1] += (cur[1]+next[1])*(cur[0]*next[1]-next[0]*cur[1]);
        }
        C[0] /= 6*A;
        C[1] /= 6*A;
        new_points.push_back(C);
    }
    return voronoi_tesselation(new_points,iter-1);
}

void save_svg(std::string filename, const std::vector<Polygon> &polygons, const std::vector<Vector> &vectors, std::string fillcol)
{
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    if (!polygons.empty()){
        for (size_t i = 0; i < polygons.size(); i++){
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \"");
            for (size_t j = 0; j < polygons[i].vertices.size(); j++){
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
    }
    if (!vectors.empty()){
        for (auto &vect : vectors){
            fprintf(f, "<g>\n");
            fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r=\"5\"/>", vect[0]*1000, 1000-vect[1]*1000);
            fprintf(f, "</g>\n");
        }
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

void save_svg(std::string filename, const std::vector<Vector> &vectors, const std::vector<Polygon> &polygons, std::string fillcol)
{
    FILE *f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (size_t i = 0; i < polygons.size(); i++){
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (size_t j = 0; j < polygons[i].vertices.size(); j++){
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    for (auto &vect : vectors){
        fprintf(f, "<g>\n");
        fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r=\"5\"/>", vect[0] * 1000, 1000 - vect[1] * 1000);
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes){
    FILE *f;
    if (frameid == 0){
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
    }
    else{
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (size_t i = 0; i < polygons.size(); i++){
        fprintf(f, "<polygon points = \"");
        for (size_t j = 0; j < polygons[i].vertices.size(); j++){
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++){
        if (frameid == j){
            fprintf(f, "inline");
        }
        else{
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++){
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1){
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}