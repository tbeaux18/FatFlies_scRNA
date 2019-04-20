#!/bin/bash


chmod a+x ./docker/build_me.sh

./docker/build_me.sh

mount_dir=$(pwd)

parentdir="$(dirname "$mount_dir")"

sudo docker run \
  --name=scrna_pipeline \
  -d \
  -v "${parentdir}":/pipeline \
  ubuntur35:pipeline
