FROM ubuntu:latest
RUN \
  apt-get update -y && \
  apt-get install -y build-essential cmake curl python-is-python3 wget\
# libssl-dev libffi-dev
    python3-dev python3 python3-pip git \
    samtools bowtie2 && \
  git clone https://github.com/Domomod/GrassSV.git && \
  pip3 install numpy pandas argcomplete quast --break-system-packages && \
  git clone https://github.com/swacisko/ALGA.git && \
  mkdir ALGA/build && \
  cd ALGA/build && \
  cmake .. && make
ENV PATH="/GrassSV:/ALGA/build:$PATH"
ENV GRASSUTILS="/GrassSV/GrassSlurm/Bash"
CMD ["/bin/bash"]