cmake_minimum_required(VERSION 3.10)
project(softrobot_project CXX)

set(CMAKE_CXX_STANDARD 17)

# Update this path to your local Eigen3 include directory if needed
include_directories(/opt/homebrew/include/eigen3)

# Source files
set(SOFTROBOT_SRC
    src/mechanics/test_softRobot.cpp
    src/mechanics/softrobot.cpp
    src/frame_util.cpp
    src/mechanics/robotState.cpp
    src/mechanics/geometry.cpp
    src/mechanics/stiffness.cpp
)

set(BEND_SRC
    src/mechanics/robotState.cpp
    src/frame_util.cpp
    src/elastics/elastic_energy.cpp
    src/elastics/bend_energy.cpp
    src/elastics/test_bend.cpp
    src/mechanics/softRobot.cpp
    src/mechanics/geometry.cpp
    src/mechanics/stiffness.cpp
)

set(STRETCH_SRC
    src/elastics/test_stretch.cpp
    src/mechanics/robotState.cpp
    src/frame_util.cpp
    src/elastics/elastic_energy.cpp
    src/elastics/stretch_energy.cpp
    src/mechanics/softRobot.cpp
    src/mechanics/geometry.cpp
    src/mechanics/stiffness.cpp
)

# Executables
add_executable(softrobot ${SOFTROBOT_SRC})
add_executable(test_bend ${BEND_SRC})
add_executable(test_stretch ${STRETCH_SRC})