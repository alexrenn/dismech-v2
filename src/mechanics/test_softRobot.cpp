#ifndef TEST_SOFTROBOT_H
#define TEST_SOFTROBOT_H

#include <iostream>
#include "softRobot.h"
#include "../eigenIncludes.h"
#include "robotState.h"
#include "environment.h"
#include "geometry.h"
#include "../frame_util.h"
#include "../springs/bend_twist_spring.h"
#include "../springs/hinge_spring.h"
#include "../springs/stretch_spring.h"
#include "../springs/triangle_spring.h"
#include "stiffness.h"

int main() {
    std::cout << "SoftRobot test initialization.\n";

    // Create GeomParams with same values
    GeomParams geom_params(1e-3,    // rod_r0
                          1e-3    // shell_h
                          );    // axs (no area scaling)

    // Create Material with same values
    Material material(1500,    // density
                     10e6,    // youngs_rod
                     10e8,    // youngs_shell
                     0.5,     // poisson_rod
                     0.30);    // poisson_shell

    // Create SimParams with same values
    SimParams dynamic_3d_sim(false,  // static_sim
                           false,    // two_d_sim
                           false,    // use_mid_edge
                           false,    // use_line_search
                           false,    // show_floor
                           true,     // log_data
                           1,        // log_step
                           1e-2,     // dt
                           25,       // max_iter
                           3.0,      // total_time
                           10,       // plot_step
                           1e-4,     // tol
                           1e-4,     // ftol
                           1e-2);    // dtol

    // Create Environment and add forces
    Environment env;
    // env.add_force("gravity", Eigen::Vector3d(0.0, 0.0, -9.81));
    // env.add_force("aerodynamics", 1.0, 10.0);  // rho=1, cd=10

    // Load geometry from file
    Geometry geo("tests/hex_parachute_n6.txt");

    // Create initial RobotState with required parameters
    Eigen::VectorXd q0 = Eigen::VectorXd::Zero(1);  // Will be resized by SoftRobot
    Eigen::VectorXd u = Eigen::VectorXd::Zero(1);   // Will be resized by SoftRobot
    std::vector<Eigen::Vector3d> a;  // Empty vector
    std::vector<Eigen::Vector3d> a1; // Empty vector
    std::vector<Eigen::Vector3d> a2; // Empty vector
    std::vector<Eigen::Vector3d> m1; // Empty vector
    std::vector<Eigen::Vector3d> m2; // Empty vector
    Eigen::VectorXi undef_ref_twist = Eigen::VectorXi::Zero(1); // Will be resized
    std::vector<Eigen::Vector3d> tau0; // Empty vector
    Eigen::VectorXi free_dof = Eigen::VectorXi::Zero(1); // Will be resized

    RobotState initial_state(q0, u, a, a1, a2, m1, m2, undef_ref_twist, tau0, free_dof);

    // Create SoftRobot instance
    SoftRobot robot(geom_params, material, geo, dynamic_3d_sim, env, initial_state);

    std::cout << "SoftRobot initialized successfully.\n";
    robot.debug();
    return 0;
}

#endif // TEST_SOFTROBOT_H