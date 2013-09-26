#!/bin/bash

source ${MODULESHOME}/init/bash

module unload PrgEnv-gnu
module unload PrgEnv-pgi
module unload PrgEnv-cray
module unload PrgEnv-intel
module unload PrgEnv-pathscale
module unload petsc/3.3-64bit-indices

module load PrgEnv-gnu
module swap gcc gcc/4.7.1

module load petsc/3.3-64bit-indices

module load acml
module load cmake
module load cray-tpsl