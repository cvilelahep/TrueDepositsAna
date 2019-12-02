#!/bin/sh                                                                                                                                                                                
DUNETPC_VERSION=v07_09_00
COMPILER=e17
DIRECTORY=dunetpc_dev_GeoEff
USERNAME=`whoami`

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export WORKDIR=${PWD}/`dirname "$0"`
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo ~`
fi

cd ${WORKDIR}
touch ${DIRECTORY}
rm -rf ${DIRECTORY}
mkdir ${DIRECTORY}
cd ${DIRECTORY}
mrb newDev -q ${COMPILER}:prof
source ${WORKDIR}/${DIRECTORY}/localProducts*/setup
mkdir work
cd srcs
mrb g -t $DUNETPC_VERSION dunetpc
cp ${WORKDIR}/TrueDepositsAnaModule  ${WORKDIR}/srcs/dunetpc/dune/TrueDepositsAna
echo "add_subdirectory(TrueDepositsAna)" >> ${WORKDIR}/srcs/dunetpc/dune/CMakeLists.txt
cd $MRB_BUILDDIR
mrbsetenv
mrb i


