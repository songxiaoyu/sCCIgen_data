
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `STsimulator` for ‘de novo’ spatial patterns

## 1. Introduction

Existing Spatially Resolved Transcriptomics (SRT) data include limited
spatial patterns. `STsimulator` generates de novel spatial patterns
using snRNAseq or single-cell SRT as reference.

`STsimulator` provides both a R package for simulating SRT data and an
interface to guide users through task selection, parameter specification
and documentation. `STsimulator` Docker

Reference: XS… ““.

## 2. Run `STsimulator` with Docker

### Installation

1.  Install `Docker Desktop`. Docker can be downloaded at
    <https://www.docker.com>.

2.  `STsimulator` uses command line prompts, such as through `Terminal`
    in Mac. You might want to download your favorite software to run
    command line. My favorite is Visual Studio Code (`VSC`). If `VSC` is
    used, Docker needs to be installed in its Extension.

3.  Clone `STsimulator` repository to your local machine, such as on
    `Terminal` type

<!-- -->

    git clone https://github.com/songxiaoyu/STsimulator

### Run Docker

1.  Open Docker Desktop.

2.  Set the path to DockerUI subdirctory

<!-- -->

    cd STsimulator/DockerUI

3.  Run the Docker container directly with a working directory (WORKDIR)
    bound to your local machine. Note WORKDIR will be the location of
    all of your input data and your outputs.

<!-- -->

    docker run --mount type=bind,source="${WORKDIR}",target=/working_directory -it songxiaoyu152/st_simulator_test

    # For Windows, one can use run.sh instead to run docker.
    sh run.sh

4.  Step 3 should open an interactive interface to guide you through
    simulation. You can provide the interface your choices of tasks,
    data, and parameters, and the `STsimulator` will provide you
    simulated data, with parameters and documentations, saved in your
    working directory.

Here are three major tasks:

#### Task 1: I want to download a pre-simulated spatial transcriptomics dataset.

`STsimulator` provides three pre-simulated datasets for researchers to
quickly check the data format and explore their analyses before the need
of spending efforts on learning the package and simulating data. These
datasets are simulated based on normal breast snRNAseq, mouse brain
SeqFISH+, and human ovarian cancer MERFISH.

Here, just select task 1 and pick the dataset to download.

#### Task 2: I want to generate a parameter file using command line prompts.

#### Task 3: I have a parameter file and want to run a simulation.

## 3. Run `STsimulator` with R package

You can install the latest version directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/RPackage/STsimulator")
library(STsimulator)
```
