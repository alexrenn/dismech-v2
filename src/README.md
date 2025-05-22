g++ -std=c++17 -I/opt/homebrew/include/eigen3 -g geometry.cpp -o geometry

g++ softrobot.cpp ../frame_util.cpp robotState.cpp -I/opt/homebrew/include/eigen3 -std=c++17 -o softrobot