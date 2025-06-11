#ifndef FRAME_UTIL_H
#define FRAME_UTIL_H

#include "eigenIncludes.h"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <cmath>

std::vector<Eigen::Vector3d> parallel_transport(
    const std::vector<Eigen::Vector3d>& u,
    const std::vector<Eigen::Vector3d>& t_start,
    const std::vector<Eigen::Vector3d>& t_end);

Eigen::Vector3d parallel_transport(const Eigen::Vector3d& u,
                                   const Eigen::Vector3d& t_start,
                                   const Eigen::Vector3d& t_end);

void compute_tfc_midedge(
    const std::vector<Eigen::Matrix3d>& p_s,                // [N faces] each 3x3
    const std::vector<Eigen::Matrix3d>& tau0_s,             // [N faces] each 3x3
    const std::vector<std::array<int, 3>>& s_s,          // [N faces] 3 signs
    std::vector<Eigen::Matrix3d>& ts,                       // output: tangents [N, 3x3]
    std::vector<std::array<double, 3>>& fs,                 // output: force projections
    std::vector<std::array<double, 3>>& cs                  // output: scalar coefficients
);

Eigen::VectorXi compute_reference_twist_util(
    const std::vector<std::array<int, 2>>& edges,
    const std::vector<std::array<int, 2>>& sgn,
    const std::vector<Eigen::Vector3d>& a1,
    const std::vector<Eigen::Vector3d>& tangent,
    const Eigen::VectorXi& ref_twist);

Eigen::Vector3d rotate_axis_angle(const Eigen::Vector3d& v,
                                    const Eigen::Vector3d& axis,
                                    double theta);

double signed_angle(const Eigen::Vector3d& u,
                    const Eigen::Vector3d& v,
                    const Eigen::Vector3d& n);

                            


#endif // FRAME_UTIL_H