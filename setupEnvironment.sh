DIRECTORY=dunetpc_dev_GeoEff
USERNAME=`whoami`

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export WORKDIR=${PWD}/`dirname "$0"`
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo ~`
fi

cd $WORKDIR/$DIRECTORY
source localProducts*/setup
# cd work
mrbslp


