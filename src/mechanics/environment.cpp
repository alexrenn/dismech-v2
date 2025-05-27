#include "environment.h"

void Environment::add_force(const Gravity& f) {
    g = f.g;
    ext_force_list.push_back("gravity");
}

void Environment::add_force(const Buoyancy& f) {
    rho = f.rho;
    ext_force_list.push_back("buoyancy");
}

void Environment::add_force(const Viscous& f) {
    eta = f.eta;
    ext_force_list.push_back("viscous");
}

void Environment::add_force(const Aerodynamics& f) {
    rho = f.rho;
    cd = f.cd;
    ext_force_list.push_back("aerodynamics");
}

void Environment::add_force(const PointForce& f) {
    pt_force = {f.x, f.y};
    pt_force_node = f.node;
    ext_force_list.push_back("pointForce");
}

void Environment::set_static() {
    if (g != 0.0f) {
        static_g = g;
    }
}

const std::vector<std::string>& Environment::get_ext_force_list() const {
    return ext_force_list;
}
