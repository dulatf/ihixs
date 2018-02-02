# iHixs 2 - Inclusive Higgs Cross Sections 
Inclusive Higgs cross sections.

## General information
iHixs is a C++ code for computing inclusive Higgs boson production cross sections in gluon fusion.

If you use this code in your publication, please cite
[arXiv:1801.XXXX](XXX).

## Dependencies
1. [LHAPDF version 6.*](https://lhapdf.hepforge.org/index.html)
2. [Boost 1.6 or higher (headers only)](http://www.boost.org/)
3. Cuba 4.2 : [a multidimensional numerical integration library](http://www.feynarts.de/cuba/)

iHixs has been succesfully built with gcc 4.8 and higher, on Linux and Mac OS X systems. With respect to C++ standards, C++11 compliance is required.

**note:** LHAPDF and Cuba should have been installed on the system with the **same compiler** as used for ihixs. This is a common source of linking problems on Mac OS, where the native compiler is CLANG LLVM.


## Download

You can obtain iHixs from the releases on the github repository,

https://github.com/dulatf/ihixs/releases

or by cloning the repository directly:

```Shell
git clone https://github.com/dulatf/ihixs.git
```

## Installation

Clone the repository or download the release tarball and unzip.
Then change create a build subdirectory in the resulting ihixs source directory and run cmake.

```Shell
cd <ihixs_src_dir>
mkdir build
cd build
cmake ..
make
```
If the dependencies cannot be found automatically run cmake instead as follows
```Shell
cmake -DLHAPDF_DIR=<lhapdf_main_dir> -DCUBA_DIR_USER=<cuba_dir> -DBOOST_DIR_USER=<boost_dir> ..
```
If you need to determine a specific set of compilers 
```Shell
      CXX=mygcc CC=mygcc cmake -DLHAPDF_DIR=<lhapdf_main_dir> -DCUBA_DIR_USER=<cuba_dir> -DBOOST_DIR_USER=<boost_dir> ..
```
## Usage


The input is controlled by a runcard. An example runcard is available
at default.card
```
    ./ihixs -i myruncard
```
Available options that can be specified in the runcard can be seen by
```
    ./ihixs --help
```
Options can also be set through command line, the gnu way, e.g.
```
    ./ihixs -i myruncard --mur=63.5
```
Some of the options have shorthands. For example 
```
    --input_filename myruncard
```
is equivalent to 
```
    -i myruncard
```  
To see for which options have available shorthands, type 
```
    ./ihixs --help 
```
and check whether an option has a shorthand after its name, for example
```    
    output_filename [o] 
```
means that the option output_filename has the shorthand 'o'



## Running the tests
* How to run tests : there are several tests integrated in the distribution, within the google test framework (https://github.com/google/googletest). You can run them by typing (from the build) 
    ```Shell
    src/tests/<name of test> 
    ```
    where <name of test> is any of the binaries in the src/test/ directory. For example 
    ```Shell
    src/tests/ihixs_eft
    ```
    runs several tests related to the higgs eft amplitudes.


## Developers
- **Achilleas Lazopoulos**
- Bernhard Mistlberger
- Falko Dulat

Made with :mortar_board: in Geneva, Menlo Park and Zurich.
