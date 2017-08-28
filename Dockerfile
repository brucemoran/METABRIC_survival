##Dockerfile to a run survival analysis on METABRIC data; 
##edit the template.metabric.R script to specify what clinical data, genes etc are required
##pass this to the Dockerfile at build with --build-arg SCRIPT=/path/to/${SCRIPT}
FROM r-base

##download
RUN R -e "install.packages('survival')"

