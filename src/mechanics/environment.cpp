#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

struct Environment {
    // Member variables for forces and properties
    float g = 0.0f;
    float rho = 0.0f;
    float eta = 0.0f;
    float cd = 0.0f;
    std::pair<float, float> pt_force = {0.0f, 0.0f};  // Example for point force (x, y)
    int pt_force_node = -1;
    float static_g = 0.0f;

    std::vector<std::string> ext_force_list;

    // Method to add forces
    void add_force(const std::string& key, const std::unordered_map<std::string, float>& kwargs) {
        if (key == "gravity") {
            g = kwargs.at("g");
        } else if (key == "buoyancy") {
            rho = kwargs.at("rho");
        } else if (key == "viscous") {
            eta = kwargs.at("eta");
        } else if (key == "aerodynamics") {
            rho = kwargs.at("rho");  // REUSING
            cd = kwargs.at("cd");
        } else if (key == "pointForce") {
            pt_force.first = kwargs.at("pt_force_x");
            pt_force.second = kwargs.at("pt_force_y");
            pt_force_node = static_cast<int>(kwargs.at("pt_force_node"));
        }
        // TODO: Translate other forces (selfContact, selfFriction, etc.)
        else {
            throw std::invalid_argument("Unknown force type");
        }

        ext_force_list.push_back(key);
    }

    // Getter for the external forces list
    const std::vector<std::string>& get_ext_force_list() const {
        return ext_force_list;
    }

    // Method to set static values
    void set_static() {
        if (g != 0.0f) {
            static_g = g;
        }
    }
};