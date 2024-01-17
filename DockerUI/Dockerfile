#FROM rpy2/base-ubuntu
#FROM python:3.9.16-slim-buster
#FROM rocker/tidyverse
FROM rocker/r-ver:4

# install R packages
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >>"${R_HOME}/etc/Rprofile.site"
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.11 python3-pip python3-setuptools python3-dev 
RUN apt-get install -y libgdal-dev
RUN apt-get install -y libfontconfig1-dev
RUN apt-get install -y libfribidi-dev
RUN apt-get install -y libharfbuzz-dev
RUN apt-get install -y libudunits2-dev

RUN pip install numpy
RUN pip install pandas
RUN pip install -U pip setuptools 
RUN apt-get install -y libbz2-dev
RUN pip install rpy2

RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('data.table', dependencies=TRUE)"
RUN R -e "install.packages('raster', dependencies=TRUE)"
RUN R -e "install.packages('spatstat', dependencies=TRUE)"
RUN R -e "install.packages('rlist', dependencies=TRUE)"
RUN R -e "install.packages('parallel', dependencies=TRUE)"
RUN R -e "install.packages('doParallel', dependencies=TRUE)"
RUN R -e "install.packages('proxy', dependencies=TRUE)"
RUN R -e "install.packages('dplyr', dependencies=TRUE)"
RUN R -e "install.packages('tibble', dependencies=TRUE)"
RUN R -e "install.packages('readr', dependencies=TRUE)"
RUN R -e "install.packages('readr', dependencies=TRUE)"
RUN R -e "install.packages('devtools', dependencies=TRUE)"
RUN R -e "devtools::install_github(\"songxiaoyu/STsimulator/RPackage\")"

# copy local files to image
RUN mkdir /working_directory
COPY docker_prompts.py docker_prompts.py
COPY py_run_r.R py_run_r.R
COPY requirements.txt requirements.txt
COPY InputData InputData
COPY RPackage RPackage
COPY ParameterFile ParameterFile
COPY example_data example_data
RUN pip install -r requirements.txt

# ENTRYPOINT 
CMD [ "python3", "docker_prompts.py" ]