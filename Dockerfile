FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    git g++ r-base r-base-dev python3 python3-pip \
    ca-certificates cpp make libltdl-dev wget unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install numpy pandas scipy 

RUN git clone --recursive https://github.com/wwylab/CliPP.git
RUN cd CliPP && python3 setup.py build