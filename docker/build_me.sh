#!/bin/bash

docker build -f ./docker/Dockerfile.ub -t ubuntu1804:base .
docker build -f ./docker/Dockerfile.r -t ubuntu:r35 .
docker build -f ./docker/Dockerfile.pipe -t ubuntur35:pipeline .
