#!/bin/bash

docker build --platform linux/amd64 -t songxiaoyu152/st_simulator_test:latest .
docker push songxiaoyu152/st_simulator_test:latest