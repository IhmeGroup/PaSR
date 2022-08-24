#include <iostream>
#include <omp.h>

#include "Solver.h"

int main(int argc, char *argv[]) {
    // Get input file from command line
    std::string input_filename = "";
    for (int i = 0; i < argc; i++){
        if (std::strncmp(argv[i],"-i",2) == 0){
            input_filename = argv[i + 1];
            std::cout << "Input file: " << input_filename << std::endl;
        }
    }
    if (input_filename.empty()){
        std::cerr << "No input file provided." << std::endl;
        throw(0);
    }

    Solver* test = new Solver(input_filename);
    test->run();
    return -1;
}