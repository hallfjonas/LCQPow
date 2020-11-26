# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Update apps on the base image
RUN apt-get -y update

# Get git, g++, make
RUN apt-get install -y git g++ make

# Clone git repo
RUN git clone https://gitlab.syscop.de/jonas.hall/lcqpOASES.git

# Change working directory
WORKDIR "lcqpOASES"

# Create directories (kept empty in git repo, thus need to create manually)
RUN mkdir bin
RUN mkdir bin/LCQP

# Compile code
RUN make

# Run example
CMD ./bin/LCQP/warm_up

