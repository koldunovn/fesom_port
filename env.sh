#!/usr/bin/env bash
#
# Source me before building/running. Mirrors fesom2/env/levante.dkrz.de/shell.gnu
# (gcc-11 + openmpi-4.1.2 + serial netcdf-c — FESOM gathers to rank 0 for I/O).
#
# Usage:
#   source /home/a/a270088/port2/fesom2_port/env.sh
# or via:
#   bash -l configure.sh

export LC_ALL=en_US.UTF-8
export CPU_MODEL=AMD_EPYC_ZEN3

module --force purge

module load git
module load gcc/11.2.0-gcc-11.2.0
module load openmpi/4.1.2-gcc-11.2.0
module load netcdf-c/4.8.1-gcc-11.2.0

export FC=mpif90 CC=mpicc CXX=mpicxx

ulimit -s unlimited
ulimit -c 0

# OpenMPI runtime knobs (verbatim from fesom2/env/levante.dkrz.de/shell.gnu).
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl=self
export OMPI_MCA_osc="pt2pt"
export UCX_IB_ADDR_TYPE=ib_global
export OMPI_MCA_coll="^ml,hcoll"
export OMPI_MCA_coll_hcoll_enable="0"
export HCOLL_ENABLE_MCAST_ALL="0"
export HCOLL_MAIN_IB=mlx5_0:1
export UCX_NET_DEVICES=mlx5_0:1
export UCX_TLS=mm,knem,cma,dc_mlx5,dc_x,self
export UCX_UNIFIED_MODE=y
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_io="romio321"
export UCX_HANDLE_ERRORS=bt
