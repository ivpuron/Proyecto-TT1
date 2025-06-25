g++ tests/test_main.cpp src/*.cpp -lm -std=c++23 -o bin/main_tests.exe
cd bin
main_tests.exe
pause