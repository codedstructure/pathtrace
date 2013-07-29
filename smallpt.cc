// adapted by Ben Bass 2013
//
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <vector>

#include "vec.h"
/*
   struct vec3 {        // Usage: time ./smallpt 5000 && xv image.ppm
   union { double x, s, r; };
   union { double y, t, g; };
   union { double z, u, b; };
   vec3(double x_=0, double y_=0, double z_=0) { x=x_; y=y_; z=z_; }
   vec3 operator+(const vec3 &b) const { return vec3(x+b.x,y+b.y,z+b.z); }
   vec3 operator-(const vec3 &b) const { return vec3(x-b.x,y-b.y,z-b.z); }
   vec3 operator*(double b) const { return vec3(x*b,y*b,z*b); }
   vec3 mult(const vec3 &b) const { return vec3(x*b.x,y*b.y,z*b.z); }
   vec3& norm() { return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
   double dot(const vec3 &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
   vec3 operator%(vec3&b) {return vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
   };*/

struct Ray { vec3 o, d; Ray(vec3 o_, vec3 d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()


struct Shape {
    virtual double intersect(const Ray& r) const = 0;
    virtual vec3 getNormal(const vec3& poi) const = 0;
    virtual Ray getLightSample(const vec3& origin, double* omega, unsigned short* Xi) const = 0;

    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    vec3 e, c;        // emission, color
protected:
    Shape(vec3 e_, vec3 c_, Refl_t refl_) :
        e(e_), c(c_), refl(refl_) {}
};

struct Sphere : public Shape {
    double rad;       // radius
    vec3 p;           // position

    Sphere(double rad_, vec3 p_, vec3 e_, vec3 c_, Refl_t refl_):
        Shape(e_, c_, refl_), rad(rad_), p(p_) {}

    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        vec3 op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps=1e-4, b=dot(op, r.d), det=b*b-dot(op, op)+rad*rad;
        if (det<0) return 0; else det=sqrt(det);
        return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
    }

    vec3 getNormal(const vec3& poi) const {
        return normalize(poi - p);
    };

    Ray getLightSample(const vec3& origin, double* omega, unsigned short* Xi) const {
        vec3 sw=p - origin;
        vec3 su=normalize((fabs(sw.x)>.1 ? vec3(0,1) : vec3(1)) % sw);
        vec3 sv=sw % su;
        double cos_a_max = sqrt(1-rad*rad/dot(origin-p, origin-p));
        double eps1 = erand48(Xi), eps2 = erand48(Xi);
        double cos_a = 1-eps1+eps1*cos_a_max;
        double sin_a = sqrt(1-cos_a*cos_a);
        double phi = 2*M_PI*eps2;
        *omega = 2*M_PI*(1-cos_a_max);
        return Ray(origin, normalize(su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a));
    };
};

// Scene: radius, position, emission, color, material
std::vector <Shape*> scene;


void buildScene() {
    scene.push_back(new Sphere(1e5, vec3( 1e5+1,40.8,81.6), vec3(),vec3(.75,.25,.25),DIFF));//Left
    scene.push_back(new Sphere(1e5, vec3(-1e5+99,40.8,81.6),vec3(),vec3(.25,.25,.75),DIFF));//Rght
    scene.push_back(new Sphere(1e5, vec3(50,40.8, 1e5),     vec3(),vec3(.75,.75,.75),DIFF));//Back
    scene.push_back(new Sphere(1e5, vec3(50,40.8,-1e5+170), vec3(),vec3(),           DIFF));//Frnt
    scene.push_back(new Sphere(1e5, vec3(50, 1e5, 81.6),    vec3(),vec3(.75,.75,.75),DIFF));//Botm
    scene.push_back(new Sphere(1e5, vec3(50,-1e5+81.6,81.6),vec3(),vec3(.75,.75,.75),DIFF));//Top
    scene.push_back(new Sphere(16.5,vec3(27,16.5,47),       vec3(),vec3(1,1,1)*.999, SPEC));//Mirr
    scene.push_back(new Sphere(16.5,vec3(73,16.5,78),       vec3(),vec3(1,1,1)*.999, REFR));//Glas
    scene.push_back(new Sphere(1.5, vec3(50,81.6-16.5,81.6),vec3(4,4,4)*100,  vec3(), DIFF)); //Lite
        /*
           Sphere(1e5, vec3(50, 1e5, 81.6), vec3(), vec3(.75, .75, .75), DIFF),
           Sphere (20, vec3(40, 20, 50), vec3(), vec3(1,0,0), DIFF),
           Sphere (20, vec3(40, 60, 50), vec3(1,1,1), vec3(), DIFF),
           */
}

inline double clamp(double x) { return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x) { return int(pow(clamp(x),1/2.2)*255+.5); }

bool intersect(const Ray &r, double &t, int &id) {
    double d;
    double inf=t=1e20;
    for(int i=scene.size(); i--;) {
        if ((d=scene[i]->intersect(r)) && d < t) {
            t=d; id=i;
        }
    }
    return t<inf;
}

vec3 radiance(const Ray &r, int depth, unsigned short *Xi, int E=1) {
    double t;                               // distance to intersection
    int id=0;                               // id of intersected object

    if (!intersect(r, t, id)) {
        return vec3(); // if miss, return black
    }

    const Shape* obj = scene[id];        // the hit object
    vec3 x=r.o+r.d*t;
    vec3 n = obj->getNormal(x);
    vec3 nl=dot(n, r.d)<0?n:n*-1;
    vec3 f=obj->c;
    double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl

    if (++depth>5 || !p) {
        if (erand48(Xi)<p) {
            f=f*(1/p);
        } else {
            return obj->e * E; //R.R
        }
    }

    if (obj->refl == DIFF) {                  // Ideal DIFFUSE reflection
        double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
        vec3 w=nl;
        vec3 u=normalize((fabs(w.x)>.1?vec3(0,1):vec3(1))%w);
        vec3 v=w%u;
        vec3 d = normalize(u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2));

        // Loop over any lights
        vec3 e;
        for (int i=0; i<scene.size(); i++){
            const Shape *s = scene[i];
            if (s->e.x<=0 && s->e.y<=0 && s->e.z<=0) continue; // skip non-lights

            double omega;
            Ray lsr = s->getLightSample(x, &omega, Xi);
            if (intersect(lsr, t, id) && id==i){  // shadow ray
                e = e + f * (s->e*dot(lsr.d, nl)*omega) * M_1_PI;  // 1/pi for brdf
            }
        }

        return obj->e*E+e+f*radiance(Ray(x,d),depth,Xi,0);
    } else if (obj->refl == SPEC)            // Ideal SPECULAR reflection
        return obj->e + f * radiance(Ray(x,r.d-n*2*dot(n, r.d)),depth,Xi);
    Ray reflRay(x, r.d-n*2*dot(n, r.d));     // Ideal dielectric REFRACTION
    bool into = dot(n, nl)>0;                // Ray from outside going in?
    double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=dot(r.d, nl), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj->e + f * radiance(reflRay,depth,Xi);
    vec3 tdir = normalize(r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:dot(tdir, n));
    double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj->e + f * (depth>2 ? (erand48(Xi)<P ?   // Russian roulette
                radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
            radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

extern "C" void addSphere(const double radius, vec3 center) {

}

extern "C" void render(const char* const fn, int w, int h, int samps) {

    buildScene();

    Ray cam(vec3(50,52,295.6), normalize(vec3(0,-0.042612,-1))); // cam pos, dir

    vec3 cx=vec3(w*.5135/h);
    vec3 cy=normalize(cx%cam.d)*.5135;
    vec3 r;
    vec3 *c=new vec3[w*h];
//#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y=0; y<h; y++) {                       // Loop over image rows
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
            for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
                for (int sx=0; sx<2; sx++) {        // 2x2 subpixel cols
                    r  = vec3();
                    for (int s=0; s<samps; s++) {
                        double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);

                        vec3 d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                            cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o+d*140,normalize(d)),0,Xi)*(1./samps);

                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + vec3(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
    }

    FILE *f = fopen(fn, "wb");         // Write image to PPM file.
    fprintf(f, "P6\n%d %d\n%d\n", w, h, 255);
    for (int i=0; i<w*h; i++) {
        fprintf(f, "%c%c%c", toInt(c[i].r), toInt(c[i].g), toInt(c[i].b));
    }
    fprintf(stderr, "\n");
}

