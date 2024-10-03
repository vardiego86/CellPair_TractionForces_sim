#!/bin/bash
# Script to create docker based runtime environment to reproduce Mpacts based simulations
# Version = 1.0
# Maintainer = Naeem Muhammad <naeem.muhammad@kuleuven.be>

DOCKER_IMAGE=mpacts_cellpir_reproducible_image
DOCKER_CONTAINER=mpacts_cellpir_reproducible_container

BASE_DIR="$PWD"		#Path to mpact-docker-reproduce-cellpair directory

docker container rm -f ${DOCKER_CONTAINER}

docker build -t ${DOCKER_IMAGE} .
docker run -di --name ${DOCKER_CONTAINER} -v "${BASE_DIR}/shared":/home/docker/data ${DOCKER_IMAGE}

docker exec -it -w /home/docker/data/cell_pair/  ${DOCKER_CONTAINER} bash