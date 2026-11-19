# Container for running the app on a server (NOT TESTED yet: no Docker where it was written).
#   docker build -t phallusia-atlas .
#   docker run -p 3838:3838 -v /path/to/phallusia.sqlite:/app/data/phallusia.sqlite:ro phallusia-atlas
# The R packages are installed in the versions recorded in app/renv.lock (the versions the app was tested with).
FROM rocker/r-ver:4.4.2
RUN R -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
COPY app/renv.lock /tmp/renv.lock
RUN R -e 'renv::restore(lockfile = "/tmp/renv.lock", library = .libPaths()[1], prompt = FALSE)'
COPY app /app
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
