FROM r-base:4.1.2

# Install the C/C++ libraries needed to run the script
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libssl-dev \
    libfontconfig1-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    build-essential \
    cmake \
    libgsl-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgeos-dev \
    libhdf5-dev \
    libboost-all-dev \
    openjdk-8-jdk \
    wget \
    git \
    libfftw3-dev \
    python3-dev python3-numpy python3-setuptools python3-matplotlib \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install the R libraries needed to run the scripts
RUN /usr/bin/R --vanilla -e 'install.packages(c("optparse", "devtools"),repos="http://cran.rstudio.com/")'
RUN /usr/bin/R --vanilla -e 'install.packages("BiocManager",repos="http://cran.rstudio.com/")'
RUN /usr/bin/R --vanilla -e 'BiocManager::install()'
RUN /usr/bin/R --vanilla -e 'BiocManager::install("Seurat")'
RUN /usr/bin/R --vanilla -e 'BiocManager::install("Signac")'
RUN /usr/bin/R --vanilla -e 'BiocManager::install(c("MOFA2", \
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",\
                       "ChIPseeker", "fastcluster", "ComplexHeatmap", "doMC"))'
RUN /usr/bin/R --vanilla -e 'BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86"))'
RUN /usr/bin/R --vanilla -e 'devtools::install_github("MarioniLab/geneBasisR")'
RUN /usr/bin/R --vanilla -e 'devtools::install_github("aertslab/cisTopic")'
RUN /usr/bin/R --vanilla -e 'install.packages("parallelDist",repos="http://cran.rstudio.com/")'
RUN /usr/bin/R --vanilla -e 'BiocManager::install("doRNG")'
RUN /usr/bin/R --vanilla -e 'BiocManager::install(c("biovizBase", "enrichTF"))'
RUN wget https://github.com/nboley/idr/archive/2.0.3.zip; \
                              unzip 2.0.3.zip; \
                              cd idr-2.0.3/; \
                              python3 setup.py install
                              
RUN /usr/bin/R --vanilla -e 'BiocManager::install(c("dorothea", "tidyverse"))' 
RUN /usr/bin/R --vanilla -e 'install.packages("pROC",repos="http://cran.rstudio.com/")' 
RUN /usr/bin/R --vanilla -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")'
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    curl

COPY . /usr/local/src/myscripts

