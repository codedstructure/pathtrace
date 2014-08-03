// adapted by Ben Bass 2013
//
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <vector>

#include <thread>
#include <future>

#include "vec.h"

struct Ray { vec3 o, d; Ray(vec3 o_, vec3 d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Rng {
    virtual double next() = 0;
};

#include "halton.h"

struct Prng : public Rng {
    Prng(long seed) {
        data[0] = 0;
        data[1] = 0;
        data[2] = seed;
    }

    double next() {
        return erand48(data);
    }
private:
    unsigned short data[3];
};

struct Shape {
    virtual double intersect(const Ray& r) const = 0;
    virtual vec3 getNormal(const vec3& poi) const = 0;
    virtual Ray getLightSample(const vec3& origin, double* omega, Rng* Xi) const = 0;

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

    Ray getLightSample(const vec3& origin, double* omega, Rng* Xi) const {
        vec3 sw=p - origin;
        vec3 su=normalize((fabs(sw.x)>.1 ? vec3(0,1) : vec3(1)) % sw);
        vec3 sv=sw % su;
        double cos_a_max = sqrt(1-rad*rad/dot(origin-p, origin-p));
        double eps1 = Xi->next(), eps2 = Xi->next();
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
    //scene.push_back(new Sphere(16.5,vec3(73,16.5,78),       vec3(),vec3(1,1,1)*.999, REFR));//Glas
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

vec3 radiance(const Ray &r, int depth, Rng *Xi, int E=1) {
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
        if (Xi->next()<p) {
            f=f*(1/p);
        } else {
            return obj->e * E; //R.R
        }
    }

    if (obj->refl == DIFF) {                  // Ideal DIFFUSE reflection
        double r1=2*M_PI*Xi->next(), r2=Xi->next(), r2s=sqrt(r2);
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
    // else REFR...
    Ray reflRay(x, r.d-n*2*dot(n, r.d));     // Ideal dielectric REFRACTION
    bool into = dot(n, nl)>0;                // Ray from outside going in?
    double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=dot(r.d, nl), cos2t;
    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj->e + f * radiance(reflRay,depth,Xi);
    vec3 tdir = normalize(r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t))));
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:dot(tdir, n));
    double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj->e + f * (depth>2 ? (Xi->next()<P ?   // Russian roulette
                radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
            radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

extern "C" void addSphere(const double radius, double x, double y, double z) {
    scene.push_back(new Sphere(radius, vec3(x, y, z), vec3(), vec3(0.5, 0.5, 0.5), DIFF));
}

extern "C" void addLight(const double radius, double x, double y, double z) {
    scene.push_back(new Sphere(radius, vec3(x, y, z), vec3(20, 20, 20), vec3(), DIFF));
}

std::mutex row_mutex;
int row_number [1];

int get_row_number() {
    std::lock_guard<std::mutex> guard(row_mutex);
    return (*row_number)++;
}

void calc_row(int y, int w, int h, Ray cam, int samps, vec3 cx, vec3 cy, vec3* c) {
#define QMC
#ifdef QMC
    Rng *Xi=new Halton(static_cast<long>(y*y*y), 17);
#else
    Rng *Xi=new Prng(static_cast<long>(y*y*y));
#endif
    for (int x=0; x<w; x++) {   // Loop cols
        for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++) {    // 2x2 subpixel rows
            for (int sx=0; sx<2; sx++) {        // 2x2 subpixel cols
                auto r = vec3();
                for (int s=0; s<samps; s++) {
                    double r1=2*Xi->next();
                    double dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                    double r2=2*Xi->next();
                    double dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);

                    vec3 d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                             cy*( ( (sy+.5 + dy)/2 + y)/h - .5) +
                             cam.d;
                    r = r + radiance(Ray(cam.o+d*140,normalize(d)),0,Xi)*(1./samps);

                } // Camera rays are pushed ^^^^^ forward to start in interior
                c[i] = c[i] + vec3(clamp(r.x),clamp(r.y),clamp(r.z)) / 4;
            }
        }
    }
}

extern "C" void render(const char* const fn, int w, int h, int samps) {
    std::vector<std::thread> thread_list;

    const Ray cam(vec3(50,52,295.6), normalize(vec3(0,-0.042612,-1))); // cam pos, dir

    const vec3 cx=vec3(w*.5135/h);
    const vec3 cy=normalize(cx%cam.d)*.5135;
    vec3 *const c=new vec3[w*h];

    int hw_threads = std::thread::hardware_concurrency();
    if (!hw_threads) {
        hw_threads = 1;
    }
    fprintf(stderr, "Threads: %d\n", hw_threads);

    int y;
    while ((y=get_row_number()) < h) {
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));

#ifdef SINGLE_THREADED
        calc_row(y,w,h,cam,samps,cx,cy,c);
        continue;
#endif
        thread_list.push_back(std::thread(calc_row, y, w, h, cam, samps, cx, cy, c));
        if (thread_list.size() >= 2*hw_threads) {
            for (auto i=0; i<hw_threads; i++) {
                thread_list.rbegin()->join();
                thread_list.pop_back();
            }
        }
    }
#ifndef SINGLE_THREADED
    while (thread_list.size()) {
        thread_list.rbegin()->join();
        thread_list.pop_back();
    }
#endif

    FILE *f = fopen(fn, "wb");         // Write image to PPM file.
    fprintf(f, "P6\n%d %d\n%d\n", w, h, 255);
    for (int i=0; i<w*h; i++) {
        fprintf(f, "%c%c%c", toInt(c[i].r), toInt(c[i].g), toInt(c[i].b));
    }
    fprintf(stderr, "\n");
}

