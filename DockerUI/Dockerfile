FROM rpy2/rpy2

# install R packages
# RUN R -e "install.packages('raster',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('spatstat',dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('rlist',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('scDesign2',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
# RUN R -e "install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')" 


# copy local files to image
RUN mkdir /working_directory
COPY simulated_data simulated_data
COPY expression_data expression_data
COPY docker_prompts.py docker_prompts.py
COPY scripts scripts

# ENTRYPOINT 
CMD [ "python3", "docker_prompts.py" ]