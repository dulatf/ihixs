# ihixs 2 - Inclusive Higgs Cross Sections
Inclusive Higgs cross sections.

# General information
iHixs is a C++ code for computing inclusive Higgs boson production cross sections in gluon fusion.

If you use this code in your publication, please cite
[arXiv:1801.XXXX](XXX).

## Download

You can obtain iHixs from the release on the github repository,

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
iHixs requires a copiler with C++11 support (e.g. gcc >= 4.8)

## Developers
- Falko Dulat
- Achilleas Lazopoulos
- Bernhard Mistlberger


Made with :mortar_board: in Geneva, Menlo Park and Zurich.
