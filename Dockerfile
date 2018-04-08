FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>

# todo: change guava to a snapshot release

# install stuff
RUN apt-get update && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip maven subversion gzip openjdk-8-jdk groovy wget \
     zlib1g-dev gnuplot lynx \
  && apt-get clean \
  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
  && cd /opt/bin && wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz \
  && gunzip KMC3.linux.tar.gz && tar -xvf KMC3.linux.tar \
#  && cd /opt && wget https://sourceforge.net/projects/bbmap/files/latest/download \
#  && mv download BBMap.tar && tar -xvzf BBMap.tar && rm BBMap.tar \
  && mkdir -p /opt/jars \
  && cd /opt && git clone https://github.com/google/guava.git \
  && cd guava/guava/ && mvn install \
#  && cd /opt && git clone https://github.com/marbl/canu.git && cd canu/src \
#  && make -j 2 \
  && cd /opt/jars \
  && wget http://mirrors.sorengard.com/apache//commons/math/binaries/commons-math3-3.6.1-bin.tar.gz \
  && gunzip commons-math3-3.6.1-bin.tar.gz && rm -f /opt/jars/commons-math3-3.6.1-bin.tar \
  && apt-get clean

# env vars
ENV NXF_OPTS "-Xms4G -Xmx20G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH

# kpi source
ADD *.nf /opt/kpi/
ADD input /opt/kpi/input/
ADD output /opt/kpi/output/
ADD src /opt/kpi/src/
ENV PATH /opt/bin:$PATH
ENV PATH /opt/bin:$PATH
ENV PATH /opt/kpi:$PATH
ENV PATH /opt/kpi/src:$PATH
ENV PATH /opt/bbmap:$PATH
ENV PATH /opt/canu/Linux-amd64/bin:$PATH
ENV CLASSPATH /opt/guava/guava/target/guava-HEAD-jre-SNAPSHOT.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:$CLASSPATH

CMD ["/bin/bash"]
#CMD ["/opt/kpi/kpil.nf"]
