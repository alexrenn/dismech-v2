#ifndef STIFFNESS_H
#define STIFFNESS_H

#include "params.h"
#include <vector>
#include <tuple>
#include <utility>
struct RodStiffness {
    double EA;
    double EI1;
    double EI2;
    double GJ;
};

struct ShellStiffness {
    std::vector<double> ks;
    double kb;
};


RodStiffness computeRodStiffness(const GeomParams& geom, const Material& material);
ShellStiffness computeShellStiffness(const GeomParams& geom, const Material& material,
                                     const std::vector<double>& ref_len, bool use_mid_edge);


#endif // STIFFNESS_H