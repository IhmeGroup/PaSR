# PaSR
Parallel partially-stirred reactor code

Some mixing model implementations are based on https://github.com/SuXY15/PaSR/blob/master/PaSR_models.py

## Table of Contents

- [PaSR](#pasr)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Clone the Repository](#clone-the-repository)
    - [Build with CMake](#build-with-cmake)
    - [Install](#install)
  - [Usage](#usage)
  - [Contributing](#contributing)
  - [References](#references)
  - [License](#license)

## Installation

### Prerequisites

Before you begin, ensure you have met the following requirements:

- C++ compiler (e.g., g++, clang++)
- [CMake](https://cmake.org/download/) - Follow the installation instructions for your platform.
- [Cantera](https://cantera.org/) - Installation instructions [here](https://cantera.org/install/index.html)
  - Make sure to install the Cantera C++ inteface!
- [OpenMP](https://www.openmp.org/) - Typically comes with most C++ compilers, but you may need to enable it during compilation.
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) - Installation instructions can be found in the project's [README](https://github.com/jbeder/yaml-cpp#building-the-code)

### Clone the Repository

Clone the project repository to your local machine:

```bash
git clone https://github.com/IhmeGroup/PaSR.git
```

### Build with CMake

Navigate to the project directory:
```bash
cd PaSR
```

Create a build directory and navigate into it:
```bash
mkdir build
cd build
```

Run CMake to configure the build:
```bash
cmake ../ -DCANTERA_ROOT=/path/to/cantera -DCMAKE_INSTALL_PREFIX=/path/to/install/bin
```

Build the project using make (or another build tool generated by CMake):
```bash
make
```

### Install

Install the compiled executable to the specified installation path:
```bash
make install
```

The executable will be installed to `/path/to/install/bin`. Optionally, you can add this directory to your PATH environment variable:
```bash
export PATH=$PATH:/path/to/install/bin
```

Now you can run the application from any directory using:
```bash
pasr -i input.toml
```

## Usage

To run the code, use the following command:
```bash
pasr -i input.toml
```
Here, `-i` specifies the input file and `input.toml` is the name of your input file in TOML format. Adjust the file name and path as needed for your specific input. Sample input files are provided in [`sample`](sample/)

## Contributing

If you'd like to contribute, please follow these steps:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature`)
3. Make your changes
4. Commit your changes (`git commit -am 'Add some feature'`)
5. Push to the branch (`git push origin feature/your-feature`)
6. Create a new pull request

## References

1. Z. Ren and S. B. Pope, “An investigation of the performance of turbulent mixing models,” Combustion and Flame, vol. 136, no. 1, pp. 208–216, Jan. 2004, doi: 10.1016/j.combustflame.2003.09.014.
2. [SuXY15/PaSR](https://github.com/SuXY15/PaSR) - A PaSR implementation in Python


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.