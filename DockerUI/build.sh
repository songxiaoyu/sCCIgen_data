#!/bin/bash

docker build --platform linux/amd64 -t songxiaoyu152/st_simulator_test .
docker push songxiaoyu152/st_simulator_test