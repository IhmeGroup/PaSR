#include <iostream>
#include <omp.h>

#include "PartiallyStirredReactor.h"

int main(int argc, char *argv[]) {
    // Get input file from command line
    std::string input_filename = "";
    for (int i = 0; i < argc; i++) {
        if (std::strncmp(argv[i], "-i", 2) == 0) {
            input_filename = argv[i + 1];
            std::cout << "Input file: " << input_filename << std::endl;
        }
    }
    if (input_filename.empty()) {
        throw std::runtime_error("No input file provided.");
    }

    PartiallyStirredReactor* reactor = new PartiallyStirredReactor(input_filename);
    reactor->initialize();
    reactor->run();
    return -1;
}