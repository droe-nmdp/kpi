FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>
# todo export TMPDIR=$PWD/work and export TMP=$PWD/work (check)
#   runs out of space in temp otherwise
# install stuff
RUN apt-get update && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip subversion gzip openjdk-8-jdk openjdk-8-doc groovy wget \
     zlib1g-dev gnuplot lynx maven \
	 bzip2 libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  && apt-get clean \
  && cd /opt && git clone https://github.com/droe-nmdp/kpi.git \
  && mkdir -p /opt/kpi/raw \
  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
  && cd /opt/bin && wget https://github.com/refresh-bio/KMC/releases/download/v3.1.1/KMC3.1.1.linux.tar.gz \
  && gunzip KMC3.1.1.linux.tar.gz && tar -xvf KMC3.1.1.linux.tar && rm -f KMC3.KMC3.1.1.linux.tar \
#  && cd /opt && wget https://sourceforge.net/projects/bbmap/files/latest/download \
#  && mv download BBMap.tar && tar -xvzf BBMap.tar && rm BBMap.tar \
  && mkdir -p /opt/jars \
#  && cd /opt && git clone https://github.com/marbl/canu.git && cd canu/src \
#  && make -j 2 \
  && cd /opt/jars \
  && wget https://www-us.apache.org/dist//commons/math/binaries/commons-math3-3.6.1-bin.tar.gz \
  && gunzip commons-math3-3.6.1-bin.tar.gz && tar -xvf commons-math3-3.6.1-bin.tar \
  && rm -f /opt/jars/commons-math3-3.6.1-bin.tar \
#  && cd /opt/ && wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 \
#  && tar -xjvf htslib-1.8.tar.bz2 && cd htslib-1.8 && make && make install \
#  && cd /opt/ && wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2 \
#  && tar -xjvf bcftools-1.8.tar.bz2 && cd bcftools-1.8 && make && make install \
#  && cd /opt/ && wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 \
#  && bunzip2 samtools-1.8.tar.bz2 \
#  && tar -xvf samtools-1.8.tar && cd samtools-1.8 && make && make install \
#  && rm -f /opt/*.bz2 /opt/*.tar* \
  && apt-get clean

# google guava
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN cd /opt && git clone https://github.com/google/guava.git \
  && cd guava/guava/ && mvn install

# env vars
ENV NXF_OPTS "-Xms4G -Xmx20G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH
ENV TMPDIR=/opt/kpi/work

# kpi source
RUN mkdir -p /opt/kpi/input /opt/kpi/output /opt/kpi/src
ENV PATH /opt/bin:$PATH
ENV PATH /opt/kpi:$PATH
ENV PATH /opt/kpi/src:$PATH
ENV PATH /opt/bbmap:$PATH
ENV PATH /opt/canu/Linux-amd64/bin:$PATH
ENV CLASSPATH /opt/guava/guava/target/guava-HEAD-jre-SNAPSHOT.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:$CLASSPATH

#CMD ["/bin/bash"]
CMD ["/opt/kpi/main.nf"]
