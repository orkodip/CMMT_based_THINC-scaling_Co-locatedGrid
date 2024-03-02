# CMMT_based_THINC-scaling_Co-locatedGrid
This repository contain the codes for CMMT-based THINC-scaling in-house CmFD solver. FVM is used to discretize the governing equations on a Co-located Grid.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ driver.cpp -o output
g++ driver_post.cpp -o post
```
To run this code

```bash
./output
./post
```
