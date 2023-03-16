#!/bin/bash

docker build --platform linux/amd64 -t annapamela/st_simulator_test .
docker push annapamela/st_simulator_test