FROM continuumio/miniconda3:4.9.2-alpine  
# tells the container to inherit from a miniconda container

WORKDIR /opt/app/   
# set the work directory

COPY . ./    
#  copy all the files and directories into the container

RUN conda env update -f environment.yml && \
    conda clean --all --yes      
# install the packages in environment.yml into containers

ENV PATH="/root/.local/bin:$PATH"      
# set the path

#RUN pip install .

RUN /opt/app/dependencies/mgl/install.sh
# command to install MGL Tools

ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"
# add MGL Tools and Autodock Vina to your path

RUN /opt/app/dependencies/adv/bin/vina --help
# Run autodock vina to ensure its properly installed
