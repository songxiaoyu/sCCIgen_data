FROM rocker/tidyverse

RUN apt-get update && apt-get install -y vim 

# Install python
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update 
RUN apt-get install -y python3.9
RUN apt-get install -y python3-pip

# # install R packages
# RUN R -e "install.packages('raster',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('spatstat',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# # RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('rlist',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('scDesign2',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')" 


# # copy local files to image
# COPY TryoutForAnna.R /bin/TryoutForAnna.R
# COPY TryoutForAnna.R /bin/TryoutForAnna.R
RUN mkdir /workdir
COPY simulated_data simulated_data
COPY expression_data expression_data
COPY docker_prompts.py docker_prompts.py

# ENTRYPOINT [ "/bin/bash", "-l", "-c" ]
CMD [ "python3", "docker_prompts.py" ]