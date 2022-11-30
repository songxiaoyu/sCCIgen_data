#!/bin/bash
BASE_DOCKER_IMAGE="st_simulator_test"

# get absolute path to working directory
read -p "Enter the absolute path to your working directory. Outputs will be saved to this directory, and all inputs must come from this directory: " -e WORKDIR

if [ "${WORKDIR:0:1}" != "/" ]
then
    echo "Please use an absolute path. It will start with '/' "
    exit
fi

docker run --mount type=bind,source="${WORKDIR}",target=/workdir -it "${BASE_DOCKER_IMAGE}"

# from within docker image, write prompts