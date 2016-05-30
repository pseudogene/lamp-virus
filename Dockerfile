#
# This file is part of lamp-virus. It provides the installation rules and commands to
# build a Docker image
#
# Copyright (C) 2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
FROM ubuntu:14.04
MAINTAINER Michael Bekaert <michael.bekaert@stir.ac.uk>

LABEL description="LAMP-VIRUS Docker" version="1.0" Vendor="Institute of Aquaculture, University of Stirling"
ENV DOCKERVERSION 1.0

USER root

COPY virus.R /opt/virus.R
COPY class_sequences.pl /usr/local/bin/class_sequences.pl
COPY collect_genomevirus.pl /usr/local/bin/collect_genomevirus.pl
COPY map_lamp.pl /usr/local/bin/map_lamp.pl

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list && \
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | apt-key add -
RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y git wget unzip
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y r-base tk tk-dev gfortran
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y libblas-dev liblapack-dev bioperl libgd-dev libgd-perl

RUN perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit'
RUN perl -MCPAN -e 'install Bio::Graphics'

RUN wget https://sourceforge.net/projects/primer3/files/primer3/2.3.7/primer3-2.3.7.tar.gz -O primer3-2.3.7.tar.gz  && \
    tar xfz primer3-2.3.7.tar.gz && \
    sed -i -e 's,/opt/primer3_config,/etc/primer3_config,g' primer3-2.3.7/src/primer3_boulder_main.c && \
    make -C primer3-2.3.7/src && \
    install primer3-2.3.7/src/primer3_core /usr/bin/ && \
    install primer3-2.3.7/src/ntdpal /usr/bin/ && \
    install primer3-2.3.7/src/ntthal /usr/bin/ && \
    install primer3-2.3.7/src/oligotm /usr/bin/ && \
    install primer3-2.3.7/src/long_seq_tm_test /usr/bin/ && \
    cp -r primer3-2.3.7/src/primer3_config /etc/primer3_config && \
    rm -rf primer3-2.3.7.tar.gz primer3-2.3.7

RUN wget http://bioinfo.unl.edu/downloads/GramAlign3_00.zip -O GramAlign3_00.zip && \
    unzip GramAlign3_00.zip && \
    make -C GramAlign3_00/src -k clean && \
    sed -i 's/CFLAGS = -O3 -funroll/CFLAGS = -O2 -funroll/' GramAlign3_00/src/Makefile && \
    make -C GramAlign3_00/src && \
    install GramAlign3_00/src/GramAlign /usr/local/bin/ && \
    rm -rf  GramAlign3_00  GramAlign3_00.zip  __MACOSX

RUN Rscript -e "install.packages('adegenet', repos='http://cran.rstudio.com', dependencies = TRUE, Ncpus = 8);"
RUN Rscript -e "install.packages('ape', repos='http://cran.rstudio.com', dependencies = TRUE, Ncpus = 8);"

RUN cd /opt && git clone https://github.com/pseudogene/lava-dna.git && \
    cd lava-dna && perl Makefile.PL && cd .. && \
    make -C lava-dna && \
    make -C lava-dna -k install

RUN mv /opt/lava-dna/t_data/loose_parameters.xml  /opt/lava-dna/t_data/looser_parameters.xml && \
    mv /opt/lava-dna/t_data/normal_parameters.xml /opt/lava-dna/t_data/loose_parameters.xml && \
    mv /opt/lava-dna/t_data/strict_parameters.xml /opt/lava-dna/t_data/normal_parameters.xml

COPY strict_parameters.xml /opt/lava-dna/t_data/strict_parameters.xml
COPY stricter_parameters.xml /opt/lava-dna/t_data/stricter_parameters.xml

RUN mkdir /virus
WORKDIR /virus
