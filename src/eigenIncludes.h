#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <memory>
#include <unistd.h>
#include <cmath>
#include <cstdlib>
#include <fcntl.h>
#include <csignal>
#include <sys/mman.h>
#include <sys/stat.h>        /* For mode constants */
#include <fcntl.h>           /* For O_* constants */
#include <numeric>           /* For computing mean and like*/

using namespace std;
using namespace Eigen;