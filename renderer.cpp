#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"



#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <omp.h>



static std::random_device rd;
static std::mt19937 rng(rd());
static std::uniform_real_distribution<double> uniform(0.,1.);


#define DIFFUSE 0
#define MIRROR 1000
#define TRANS 2000



double I = 2e5;
double eps = 1e-8;
class Vector {      
    public:
        double coords[3];
        explicit Vector(double x = 0., double y = 0., double z = 0.){
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
            };
        Vector& operator +=(const Vector& b ) {
            coords[0] += b[0];
            coords[1] += b[1];
            coords[2] += b[2];
            return *this;
            };
        
        Vector& operator *=(double a) {
            coords[0] *= a;
            coords[1] *= a;
            coords[2] *= a;
            return *this;
            };
        const double& operator[](int i) const {return coords[i];}
        double& operator[](int i) {return coords[i];}
    };

Vector operator+(const Vector& a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
    };

Vector operator-(const Vector& a, const Vector &b){
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
};

Vector operator*(double a, const Vector& V){
    return Vector(a*V[0], a*V[1], a*V[2]);
    };
Vector operator*(const Vector& V, double a){
    return Vector(a*V[0], a*V[1], a*V[2]);
    };

Vector operator*(const Vector& U, const Vector& V){
    return Vector(U[0]*V[0], U[1]*V[1], U[2]*V[2]);
    };

Vector operator/(const Vector& V, double a){
    return (1/a)*V;
};

Vector abs(const Vector &V){
    return Vector(abs(V[0]),abs(V[1]),abs(V[2]));
};





class Light{
    public:
        Vector origin;
        double intensity;
        explicit Light(Vector C = Vector(0, 0, 0), double I = pow(10, 5))
        {
            origin = C;
            intensity = I;
        }
};


double min_value(const Vector& N){
    std::vector<double> v;
    v.push_back(N[0]);
    v.push_back(N[1]);
    v.push_back(N[2]);
    return *std::min_element(v.begin(), v.end());
};
int min_index(const Vector &N){
    std::vector<double> v;
    v.push_back(N[0]);
    v.push_back(N[1]);
    v.push_back(N[2]);
    return std::min_element(v.begin(),v.end()) - v.begin();
}

double dot(const Vector& a, const Vector& b ) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    };

Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]); 
}

double norm(const Vector& a){
    return sqrt(dot(a,a));
    };

Vector normalize(const Vector& a){
    return (1/norm(a))*a;
};

Vector power(const Vector& a, double power){
    return Vector(pow(a[0],power), pow(a[1],power), pow(a[2],power));
}



typedef struct Ray{
    Vector O;
    Vector u;
    double refr;
    //Ray(Vector O, Vector u, double refr){O=O; u=u; refr=refr;};
} Ray;

typedef struct Intersection{
    bool intersect;
    Vector P;
    Vector N;
    Vector u;
    int id;
    double t;
} Intersection;


class Sphere{
    public:
        Vector C;
        double R;
        Vector albedo;
        int mat;
        double refr;
        Sphere(Vector center, double radius, Vector colour, int m, double r){
            this->C = center;
            this->R = radius;
            this->albedo = colour;
            this->mat = m;
            this->refr = r;
        };
        
        Intersection intersect(Ray &ray){
            Intersection inter;
            double t,t1,t2;
            Vector u = ray.u;
            Vector O = ray.O;
            Vector OmC = O-C;
            double D = pow(dot(u,OmC),2) - (dot(OmC,OmC) - R*R);

            if(D < 0){
                inter.intersect = false;
            }
            else{
                t1 = dot(u,C-O) - sqrt(D);
                t2 = dot(u,C-O) + sqrt(D);
                if(t2 < 0){
                    inter.intersect = false;
                }
                else{
                    inter.intersect = true;
                    if(t1 >= 0){
                        t = t1;
                    }
                    else{
                        t = t2;
                    }
                    inter.t = t;
                    inter.P = O + t*u;
                    inter.N = normalize(inter.P - C);
                    inter.u = u;
                }
            }
            return inter;
        }
};

double uni() { return uniform(rng); }
/*
void box_muller(double stdev, double& x, double& y) {
    double r = sqrt(-2*log(uni()))*stdev;
    double ang = 2*M_PI * uni();
    x = r*cos(ang);
    y = r*sin(ang);
}

*/

class Scene {
    public:
        std::vector<Sphere> Spheres;
        Vector Camera;
        Ray Source;
        double refr;
        bool fresnel;
        Light light;
        Scene(std::vector<Sphere> S, Vector c, Ray s, double r, bool fres, Light l){
            this->Camera = c;
            this->Spheres = S;
            this->Source = s;
            this->refr = r;
            this->fresnel = fres;
            this->light = l;

        };
        Intersection intersect(Ray &ray){
            Intersection sceneinter;
            double minimum = std::numeric_limits<double>::infinity();
            double d;
            Intersection inter;
            inter.intersect = false;
            for(int i = 0; i<Spheres.size(); i++){
                sceneinter = Spheres[i].intersect(ray);
                if(sceneinter.intersect and (d = norm(sceneinter.P - ray.O)) < minimum){
                    minimum = d;
                    //inter.intersect = true;
                    inter = sceneinter;
                    inter.id = i;
                    };
                };
            return inter;
        }

    Vector random_cos(const Vector& N) {
    double z = uni(), ang = 2*M_PI*uni();
    double r = sqrt(1-z);
    double x = r * cos(ang);
    double y = r * sin(ang);
    z = sqrt(z);
    int index = min_index(abs(N));
    Vector T1 = N;
    T1[index] = 0;
    switch(index){
        case 0:
        T1[1] = N[2], T1[2] = -N[1];
        break;
        case 1:
        T1[0] = N[2], T1[2] = -N[0];
        break;
        case 2:
        T1[0] = N[1], T1[1] = -N[0];
        break;
    }
    T1 = normalize(T1);
    Vector T2 = cross(N,T1);
    return x*T1+ y*T2 + z*N;
}

        Vector getcolor(Ray &ray, int ray_depth){
            if(ray_depth < 0) return Vector(0.,0.,0.);
            Intersection sceneinter = intersect(ray);
            Vector colour(0.,0.,0.);
            
            //if(sceneinter.intersect){
            Sphere s = Spheres[sceneinter.id];
            Vector newP = sceneinter.P + eps*sceneinter.N;
            Vector newP2 = sceneinter.P - eps*sceneinter.N;
            //if(s.mat != MIRROR && s.mat != TRANS && s.mat != DIFFUSE){ printf("%u, id = %d, %f %f %f\n",s.mat,sceneinter.id,s.albedo[0],s.albedo[1],s.albedo[2]);}
            switch(s.mat){

                case MIRROR:{
                    Ray reflected_ray;
                    reflected_ray.O = newP;
                    reflected_ray.u = sceneinter.u - 2*dot(sceneinter.u,sceneinter.N)*sceneinter.N;
                    reflected_ray.refr = refr;
                    colour = getcolor(reflected_ray, ray_depth-1);
                    break;}

                case TRANS:{
                    Vector N = sceneinter.N;
                    Vector P = sceneinter.P + 0.01*N;
                    P = newP;
                    double n1;
                    double n2;

                    if (dot(ray.u, N) > 0){
                        N = (-1)*N;
                        //P = sceneinter.P + 0.01*N;
                        P = newP2;
                        
                        n1 = s.refr; 
                        n2 = 1;
                        
                    }else{
                        n1 = 1;
                        n2 = s.refr;
                    }

                    double n1n2 = n1 / n2;
                    double dot_dir = dot(ray.u, N);
                    double radicand = 1 - pow(n1n2, 2) * (1 - pow(dot_dir, 2));

                    double k0 = pow(n1 - n2, 2) / pow(n1 + n2, 2);
                    double R = k0 + (1 - k0) * pow(1 - abs(dot(N,ray.u)), 5);
                
                    //double u = uniform(engine_scene[omp_get_thread_num()]);
                    double u = uni();


                    if (u>R && radicand >= 0){ 
                        Vector r_tan = n1n2 * (ray.u - dot_dir * N);
                        Vector r_nor = (-sqrt(radicand))*N;

                        Ray refract;
                        refract.O = P - 2*eps*N;
                        refract.u = r_tan + r_nor;

                        return getcolor(refract, ray_depth - 1);
                    }
                    Ray reflect;
                    reflect.O = P;
                    reflect.u = ray.u - 2 * dot_dir * N;
                    return getcolor(reflect, ray_depth - 1);
                    break;
                       
                }
                case DIFFUSE:{                    
                    Source.u = normalize(newP - Source.O);
                    Intersection inter = intersect(Source);
                    double d;
                    if(norm(newP - Source.O) <= norm(inter.P-Source.O)){
                        d = norm(Source.O - newP);
                        //colour = Spheres[sceneinter.id].albedo;
                        colour = (I*std::max(dot(inter.N, normalize(Source.O - newP)),0.)/(4*M_PI*M_PI*d*d)) * (Spheres[inter.id].albedo);
                        //colour[0] = std::min(255., pow(colour[0], 1./2.2)*255);
                        //colour[1] = std::min(255., pow(colour[1], 1./2.2)*255);
                        //colour[2] = std::min(255., pow(colour[2], 1./2.2)*255);

                    }
                    Ray randomRay;
                    randomRay.O = inter.P;
                    randomRay.u = random_cos(sceneinter.N);
                    randomRay.refr = refr;
                        //Vector col2 = 
                        //printf("%f %f %f \n", col2[0],col2[1],col2[2]);

                    colour += s.albedo * getcolor(randomRay, ray_depth-1);;
                    
    
                    break;
                    };
            }
            
        return colour;
        }
};


int main() {
    
    Vector C1(0,1000,0);
    Vector Red(1,0,0);
    Vector C2(0,0,-1000);
    Vector Green(0,1,0);
    Vector C3(0,-1000,0);
    Vector Blue(0,0,1);
    Vector C4(0,0,1000);
    Vector Purple(1,0,1);
    Vector C5(-1000,0,0);
    Vector Aqua(0,1,1);
    Vector C6(1000,0,0);
    Vector Yellow(1,1,0);    

    
    Vector Grey(0.5,0.5,0.5);

    Sphere RedSphere(C1,940,Red,DIFFUSE,1);
    Sphere GreenSphere(C2,940,Green,DIFFUSE,1);

    Sphere BlueSphere(C3,990,Blue,DIFFUSE,1);

    Sphere PurpleSphere(C4,940,Purple,DIFFUSE,1);
    Sphere AquaSphere(C5,940,Aqua,DIFFUSE,1);
    Sphere YellowSphere(C6,940,Yellow, DIFFUSE,1);

    Sphere MirrorSphere(Vector(-20,0,0), 10, Grey, MIRROR, 0);
    Sphere GreySphere(Vector(0,0,0), 10, Grey, DIFFUSE, 1.5);
    Light L = Light(Vector(-10,20,40),pow(10,5));
    


    //Sphere middle_sphere(Vector(0, 0, 0), Vector(0.5, 0.5, 0.5), 10, false, true, 1.5);
    Sphere right_sphere(Vector(20, 0, 0), 10, Vector(0.5, 0.5, 0.5), TRANS, 1.5);
    //Sphere right_hollow_sphere(Vector(20, 0, 0), Vector(0.5, 0.5, 0.5), 10 - 1, false, true, 1./1.5);

    std::vector<Sphere> Spheres1;
    
    Spheres1.push_back(YellowSphere);
    Spheres1.push_back(AquaSphere);
    Spheres1.push_back(RedSphere);
    Spheres1.push_back(GreenSphere);
    Spheres1.push_back(BlueSphere);
    Spheres1.push_back(PurpleSphere);

    Spheres1.push_back(GreySphere);
    //Spheres1.push_back(MirrorSphere);
    //Spheres1.push_back(right_sphere);
    
    
    Ray Source;
    Source.O = Vector(-10,20,40);


    int Qx = 0;
    int Qy = 0;
    int Qz = 55;

    Vector Camera(Qx,Qy,Qz);

    bool fresnel = true;
    Scene scene(Spheres1, Camera, Source,1,fresnel,L);

    int rays = 36;
    //if (fresnel) rays = 600;
    int W = 512;
    int H = 512;
    int x;
    int y;
    float alpha = M_PI/3;
    Vector direction(0,0,0);
    Ray ray;
    ray.O = Camera;
    ray.u = direction;
                    
    ray.refr = 1;
    double ray_depth = 5;
    

    std::vector<unsigned char> image(W*H * 3, 0);
    #pragma omp parallel for, schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        //printf("%d\n",i);
        for (int j = 0; j < W; j++) {
            x = j;
            y = H-i-1;

            direction[0] = Qx + x + 0.5 - W/2;
            direction[1] = Qy + y + 0.5 - H/2;
            direction[2] = Qz - W/(2*tan(alpha/2));
            Vector color = Vector(0,0,0);
            for (int i = 0; i < rays;i++){
                ray.u = normalize(direction - Camera);

                color += scene.getcolor(ray, ray_depth);
            };
            color = color/rays;
            //image[(i*W + j) * 3 + 0] = std::min(255.,pow(color[0],1./2.2));
            //image[(i*W + j) * 3 + 1] = std::min(255.,pow(color[1],1./2.2));
            //image[(i*W + j) * 3 + 2] = std::min(255.,pow(color[2],1./2.2));
            color[0] = std::min(255., pow(color[0], 1./2.2)*255);
            color[1] = std::min(255., pow(color[1], 1./2.2)*255);
            color[2] = std::min(255., pow(color[2], 1./2.2)*255);

            image[(i*W + j) * 3 + 0] = color[0];
            image[(i*W + j) * 3 + 1] = color[1];
            image[(i*W + j) * 3 + 2] = color[2];
            
        }
    
    }
    stbi_write_png("render2.png", W, H, 3, &image[0], 0);
 
    return 0;

}