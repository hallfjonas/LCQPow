##
##	This file is part of LCQPanther.
##
##	LCQPanther -- A Solver for Quadratic Programs with Commplementarity Constraints.
##	Copyright (C) 2020 - 2021 by Jonas Hall et al.
##
##	LCQPanther is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	LCQPanther is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with LCQPanther; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##


# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Update apps on the base image
RUN apt-get -y update && apt-get -y upgrade

# Create time zone variables for tzdata config
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Get git, g++, make
RUN apt-get install -y gcc g++ gfortran subversion git make cmake liblapack-dev pkg-config --install-recommends

# Get qpOASES and CasADi
RUN svn co https://projects.coin-or.org/svn/qpOASES/stable/3.2 /usr/local/qpOASES
RUN git clone --depth 1 --branch 3.5.5 https://github.com/casadi/casadi.git /usr/local/casadi

# Build qpOASES
WORKDIR "/usr/local/qpOASES"
RUN make

# Build CasADi
RUN mkdir "/usr/local/casadi/build"
WORKDIR "/usr/local/casadi/build"
RUN cmake ..
RUN make
RUN make install

# Clone git repo
RUN git clone https://gitlab.syscop.de/jonas.hall/LCQPanther.git -b stable /usr/local/LCQPanther

# Change working directory
WORKDIR "/usr/local/LCQPanther"
RUN make

# Run example
CMD bin/warm_up
