#ifndef INCLUDE_SMALLPT_H
#define INCLUDE_SMALLPT_H

extern "C" void addSphere(double radius, double x, double y, double z);
extern "C" void addLight(double radius, double x, double y, double z);
void buildScene();
extern "C" void render(const char* const fn, int w, int h, int samps);

#endif
