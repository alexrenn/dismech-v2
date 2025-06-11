#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include <iostream>
using namespace std;
#include "../eigenIncludes.h"


struct GeometryData {
  const std::vector<Eigen::Vector3d> nodes;
  const std::vector<std::array<int, 2> > edges;
  const std::vector<std::array<int, 3> > face_nodes;
};

struct GeomParams {
    double rod_r0;
    double shell_h;
    double* axs = nullptr;
    double* jxs = nullptr;   // Optional, initialized to nullptr
    double* ixs1 = nullptr;  // Optional, initialized to nullptr
    double* ixs2 = nullptr;  // Optional, initialized to nullptr

    GeomParams(double rod_r0, double shell_h, double* axs = nullptr, double* jxs = nullptr,
               double* ixs1 = nullptr, double* ixs2 = nullptr)
        : rod_r0(rod_r0), shell_h(shell_h), axs(axs), jxs(jxs), ixs1(ixs1), ixs2(ixs2) {}
};

struct Material {
    double density;
    double youngs_rod;
    double youngs_shell;
    double poisson_rod;
    double poisson_shell;

    Material(double density, double youngs_rod, double youngs_shell,
             double poisson_rod, double poisson_shell)
        : density(density), youngs_rod(youngs_rod), youngs_shell(youngs_shell),
          poisson_rod(poisson_rod), poisson_shell(poisson_shell) {}
    //std::cout << "\npoisson_shell INSIDE PARAMS: " << poisson_shell << std::endl;
};

struct SimParams {
    bool static_sim;
    bool two_d_sim;
    bool use_mid_edge;
    bool use_line_search;
    bool log_data;
    int log_step;
    bool show_floor;
    double dt;
    int max_iter;
    double total_time;
    int plot_step;
    double tol;
    double ftol;
    double dtol;
    std::string solver = "np";  // Default value
    bool sparse = false;

    SimParams(bool static_sim, bool two_d_sim, bool use_mid_edge, bool use_line_search,
              bool log_data, int log_step, bool show_floor, double dt, int max_iter,
              double total_time, int plot_step, double tol, double ftol, double dtol,
              const std::string& solver = "np", bool sparse = false)
        : static_sim(static_sim), two_d_sim(two_d_sim), use_mid_edge(use_mid_edge),
          use_line_search(use_line_search), log_data(log_data), log_step(log_step),
          show_floor(show_floor), dt(dt), max_iter(max_iter), total_time(total_time),
          plot_step(plot_step), tol(tol), ftol(ftol), dtol(dtol), solver(solver), sparse(sparse) {}
};

#endif // PARAMS_H