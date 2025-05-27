#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include <string>
#include <vector>
#include <utility>

struct Gravity {
    float g;
};

struct Buoyancy {
    float rho;
};

struct Viscous {
    float eta;
};

struct Aerodynamics {
    float rho;
    float cd;
};

struct PointForce {
    float x;
    float y;
    int node;
};

class Environment {
public:
    // Member variables
    float g = 0.0f;
    float rho = 0.0f;
    float eta = 0.0f;
    float cd = 0.0f;
    std::pair<float, float> pt_force = {0.0f, 0.0f};
    int pt_force_node = -1;
    float static_g = 0.0f;

    std::vector<std::string> ext_force_list;

    // Force addition overloads
    void add_force(const Gravity& f);
    void add_force(const Buoyancy& f);
    void add_force(const Viscous& f);
    void add_force(const Aerodynamics& f);
    void add_force(const PointForce& f);

    // Utility
    void set_static();
    const std::vector<std::string>& get_ext_force_list() const;
};

#endif // ENVIRONMENT_H