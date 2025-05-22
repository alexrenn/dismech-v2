#ifndef FRAME_UTIL_H
#define FRAME_UTIL_H

#include "eigenIncludes.h"
#include <vector>

std::vector<Eigen::Vector3d> parallel_transport(
    const std::vector<Eigen::Vector3d>& u,
    const std::vector<Eigen::Vector3d>& t_start,
    const std::vector<Eigen::Vector3d>& t_end);

Eigen::Vector3d parallel_transport(const Eigen::Vector3d& u,
                                   const Eigen::Vector3d& t_start,
                                   const Eigen::Vector3d& t_end);

#endif // FRAME_UTIL_H