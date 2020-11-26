# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Update apps on the base image
RUN apt-get -y update

# Create time zone variables for tzdata config
ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Get git, g++, make
RUN apt-get install -y gcc g++ gfortran git make cmake liblapack-dev pkg-config --install-recommends

# Get and build CasADi
RUN git clone https://github.com/casadi/casadi.git -b master casadi
RUN mkdir build
WORKDIR "/casadi/build"
RUN cmake ..
RUN make
RUN make install
RUN echo "." >> /etc/ld.so.conf
RUN ldconfig
WORKDIR "/"

# Clone git repo
RUN git clone https://gitlab.syscop.de/jonas.hall/lcqpOASES.git

# Change working directory
WORKDIR "/lcqpOASES"

# Create directories (kept empty in git repo, thus need to create manually)
RUN mkdir bin
RUN mkdir bin/LCQP

# Compile code
RUN make

# Run example
CMD ./bin/LCQP/warm_up

