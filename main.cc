#include "stdlib.h"
#include "smallpt.h"


int main(int argc, char *argv[]) {
    int w=512, h=384, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
    render("image.ppm", w, h, samps);
}
