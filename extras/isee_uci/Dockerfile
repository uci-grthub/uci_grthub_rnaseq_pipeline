# This docker image has LaTeX to build the vignettes
FROM bioconductor/shiny

# Update apt-get
RUN apt-get update && apt-get upgrade -y && apt-get clean && rm -rf /var/lib/apt/lists/*    

RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils \
	&& apt-get install -y --no-install-recommends \
	vim \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN R -e 'BiocManager::install("BiocStyle")'

## Install BiocStyle
RUN R -e 'BiocManager::install("iSEE")'
RUN R -e 'BiocManager::install("iSEEde")'
RUN R -e 'BiocManager::install("airway")'
RUN R -e 'BiocManager::install("DESeq2")'
