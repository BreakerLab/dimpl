# Copyright (c) 2018 Kenneth Brewer.

FROM jupyter/scipy-notebook:abdb27a6dfbb

LABEL maintainer="Kenneth Brewer <kenibrewer@gmail.com>"

# Install any necessary packages with apt-get
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxtst6 \
    libgconf2-4 \
    libxss1 \
    gnuplot-x11 \
    cufflinks \
    zlib1g-dev \
    vim \
    less \
    xvfb && \
    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Download orca AppImage, extract it, and make it executable under xvfb
RUN wget https://github.com/plotly/orca/releases/download/v1.1.1/orca-1.1.1-x86_64.AppImage -P /home
RUN chmod 777 /home/orca-1.1.1-x86_64.AppImage

# To avoid the need for FUSE, extract the AppImage into a directory (name squashfs-root by default)
RUN cd /home && /home/orca-1.1.1-x86_64.AppImage --appimage-extract
RUN printf '#!/bin/bash \nxvfb-run --auto-servernum --server-args "-screen 0 640x480x24" /home/squashfs-root/app/orca "$@"' > /usr/bin/orca
RUN chmod 777 /usr/bin/orca
RUN chmod -R 777 /home/squashfs-root/

# Install RNAcode
RUN cd /opt && \
    wget http://github.com/downloads/wash/rnacode/RNAcode-0.3.tar.gz && \
    tar -xzf RNAcode-0.3.tar.gz && \
    rm RNAcode-0.3.tar.gz && \
    cd RNAcode-0.3 && \
    ./configure && \
    make && \
    make install

RUN cd /opt && \
    wget http://eddylab.org/software/rscape/rscape.tar.gz && \
    tar -xzf rscape.tar.gz && \
    rm rscape.tar.gz &&\
    cd rscape_v* && \
    ./configure && \ 
    make && \
    make install && \ 
    cp bin/* /usr/local/bin

#Install Turbofold
RUN cd /opt/ && wget http://rna.urmc.rochester.edu/Releases/current/RNAstructureLinuxTextInterfaces64bit.tgz && \
    tar xf RNAstructureLinuxTextInterfaces64bit.tgz && rm RNAstructureLinuxTextInterfaces64bit.tgz
RUN ln -fs /opt/RNAstructure/exe/* /usr/local/bin/

# Install any necessary conda packages
USER ${NB_USER}
RUN rm /opt/conda/conda-meta/pinned
RUN conda install --quiet --yes \
    'ipywidgets>7.5' \
    'conda-forge::python-dotenv' \
    'biopython' \
    'bioconda::bedops' \
    'bioconda::bx-python' \
    'bioconda::bedtools' \
    'bioconda::pybedtools' \
    'bioconda::pysam' \
    'plotly::plotly' \
    'plotly::plotly-orca' \
    'mysql-connector-python'  && \
    conda clean --all -f -y

# Install Genomeview and add fix for Defective Importerrors
RUN pip install genomeview==1.0.1
ADD genomeview-1.0.1_dimpl.patch /opt/conda/lib/python3.7/site-packages/genomeview/
RUN cd /opt/conda/lib/python3.7/site-packages/genomeview/ && \
    git apply genomeview-1.0.1_dimpl.patch && \
    rm genomeview-1.0.1_dimpl.patch

# Install Globus-cli
RUN pip install globus-cli

USER root
RUN fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER
