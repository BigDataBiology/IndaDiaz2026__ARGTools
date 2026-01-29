# Use a multi-architecture base image
FROM rocker/shiny:4.4

ENV DEBIAN_FRONTEND=noninteractive

# Install essential system dependencies and R
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages

#RUN R -e "install.packages(c('bslib', 'dplyr', 'ggplot2', 'gridExtra', 'tidyverse', 'RColorBrewer', 'reactable'), repos='https://cloud.r-project.org/')"

RUN R -e "install.packages('pak')"
RUN R -e "options(repos = c(CRAN='https://packagemanager.posit.co/cran/latest')); pak::pkg_install(c('bslib','dplyr','ggplot2','gridExtra','RColorBrewer','reactable'))"

#RUN R -e "options(repos = c(CRAN='https://packagemanager.posit.co/cran/latest')); #pak::pkg_install(c('tidyverse'))"


# Copy your Shiny app into the container
COPY app/ /srv/shiny-server/argcompare_app/

# Set permissions so Shiny Server can access the app
RUN chown -R shiny:shiny /srv/shiny-server && chmod -R 755 /srv/shiny-server
RUN rm -r /srv/shiny-server/0*
RUN rm -r /srv/shiny-server/1*
RUN rm -r /srv/shiny-server/index.html


# Expose Shiny Server port
EXPOSE 3838

# Start Shiny Server
CMD ["/usr/bin/shiny-server"]
