#!/bin/bash
BASE_DOCKER_IMAGE="songxiaoyu152/st_simulator_test"

printf "\t-----------Start ST simulator-----------\n"
# get absolute path to working directory
read -p "Enter the absolute path to your working directory. Outputs will be saved to this directory, and all inputs must come from this directory: " -e WORKDIR

if [ "${WORKDIR:0:1}" != "/" ]
then
    echo "Please use an absolute path. It will start with '/' "
    exit
fi

docker run --mount type=bind,source="${WORKDIR}",target=/working_directory -it "${BASE_DOCKER_IMAGE}"

# from within docker image, write prompts