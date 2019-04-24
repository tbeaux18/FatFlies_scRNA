#!/bin/bash


mount_dir=$(pwd)

chmod a+x "${mount_dir}"/docker/build_me.sh

./docker/build_me.sh

parentdir="$(dirname "$mount_dir")"

sudo docker run \
  --name=scrna_pipeline \
  -d \
  -v "${parentdir}":/pipeline \
  ubuntur35:pipeline


# mv the necessary output files one directory up
# grab relevant mapping results
# check the yaml log file to make sure output code is [1] 0 indicates nothing is wrong
# once the genome index is built, move one directory up so resources are not wasted on redoing it
