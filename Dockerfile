FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>
# todo: change the supercsv stuff

# install stuff; mainly
# https://github.com/medvedevgroup/bloomtree-allsome
RUN apt-get update && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip maven subversion gzip openjdk-8-jdk groovy wget \
     zlib1g-dev \
  && apt-get clean \
  && cd /opt && mvn dependency:get -Dartifact=net.sf.supercsv:super-csv:2.4.0 \
  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
  && cd /opt/ && git clone https://github.com/medvedevgroup/bloomtree-allsome \
  && cd bloomtree-allsome && export HOME=/opt \
  && mkdir deps && cd deps && set -e \
  && wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz \
  && tar xf jellyfish-2.2.6.tar.gz && rm jellyfish-2.2.6.tar.gz && cd jellyfish-2.2.6 \
  && ./configure --prefix=/opt && make -j 4 && make install && cd .. \
  && git clone --recursive https://github.com/simongog/sdsl-lite.git \
  && cd sdsl-lite && ./install.sh && cd .. \
  && git clone https://github.com/RoaringBitmap/CRoaring.git \
  && cd CRoaring && mkdir build && cd build \
  && cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt ..  \
  && make install && export LD_LIBRARY_PATH=$HOME/lib:\$LD_LIBRARY_PATH \
  && cd /opt/bloomtree-allsome/src && make && cd ../bfcluster && make \
  && apt-get clean

# env vars
ENV NXF_OPTS "-Xms1G -Xmx5G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH

# kharsh source
ADD *.nf /opt/kpi/
ADD input /opt/kpi/input/
ADD src /opt/kpi/src/
ENV PATH /opt/bin:$PATH
ENV PATH /opt/kpi:$PATH
ENV PATH /opt/kpi/src:$PATH
ENV PATH bloomtree-allsome/src/bfcluster:$PATH
ENV CLASSPATH /opt/kpi/src/super-csv.jar:$CLASSPATH

CMD ["/bin/bash"]
#CMD ["/opt/kpi/kpil.nf"]