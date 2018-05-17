# Use an official Python runtime as a parent image
FROM ubuntu:14.04

# Generate and set locale to support unicode
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV PYTHONIOENCODING UTF-8

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
ADD DEPENDENCIES /app

RUN apt-get -y -q update
RUN apt-get -y -q install wget
RUN apt-get -y -q install gcc  # We need to install C libs to install PSI4/gcc-5 dependency
RUN apt-get -y -q install libxext6 libxrender1  # Necessary for smallmol tests (rdkit)

# Install Miniconda
ENV CONDA_INSTALLER Miniconda3-latest-Linux-x86_64.sh
ENV CONDA_PATH /home/user/conda
RUN cd ~ &&\
    wget https://repo.continuum.io/miniconda/$CONDA_INSTALLER 2> /dev/null &&\
    bash $CONDA_INSTALLER -b -p $CONDA_PATH &&\
    rm $CONDA_INSTALLER
ENV PATH $CONDA_PATH/bin:$PATH
ENV LD_LIBRARY_PATH $CONDA_PATH/lib

RUN conda config --set always_yes yes
RUN conda update -n base conda

RUN conda install conda conda-build anaconda anaconda-client requests pytest -y -q  # this has to go into the root environment

RUN conda create -q -n travis-env python=$TRAVIS_PYTHON_VERSION
RUN /bin/bash -c "source activate travis-env"

RUN conda config --add channels acellera
RUN conda config --add channels psi4
RUN conda install anaconda-client -y -q
RUN conda update --all -y -q

# Install the Conda dependencies HTMD needs
RUN conda install --file DEPENDENCIES -y

# Cleaning to reduce amount of space used on disk
RUN conda clean --all -y 

# CMD /bin/bash -c "source activate travis-env"
CMD bash