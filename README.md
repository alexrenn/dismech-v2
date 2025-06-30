# The following tests require manual comparison with dismech python.

# Method 1: With CMake 
Navigate to /build. Use commands to compile and build.

cmake ..
make

# To run, use the following commands
./softrobot
./test_bend
./test_stretch

# Method 2: Manual (Without CMake)

# To test softRobot and Geometry class, use test_softRobot.cpp and load appropriate input file. Replace -I/opt/homebrew/include/eigen3 to local location of library.

# To compile
g++ src/mechanics/test_softRobot.cpp src/mechanics/softrobot.cpp src/frame_util.cpp src/mechanics/robotState.cpp src/mechanics/geometry.cpp src/mechanics/stiffness.cpp -I/opt/homebrew/include/eigen3 -std=c++17 -o softrobot

# To run
./softRobot

# To test bending class
g++ src/mechanics/robotState.cpp src/frame_util.cpp src/elastics/elastic_energy.cpp src/elastics/bend_energy.cpp src/elastics/test_bend.cpp src/mechanics/softRobot.cpp src/mechanics/geometry.cpp src/mechanics/stiffness.cpp -I/opt/homebrew/include/eigen3 -std=c++17 -o test_bend

# To test stretching class
g++ src/elastics/test_stretch.cpp src/mechanics/robotState.cpp src/frame_util.cpp src/elastics/elastic_energy.cpp src/elastics/stretch_energy.cpp src/mechanics/softRobot.cpp src/mechanics/geometry.cpp src/mechanics/stiffness.cpp -I/opt/homebrew/include/eigen3 -std=c++17 -o test_stretch
