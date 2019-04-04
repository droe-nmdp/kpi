FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>

# apt stuff
RUN apt-get update \
  && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip subversion gzip openjdk-8-jdk openjdk-8-doc groovy wget \
     zlib1g-dev gnuplot lynx maven \
     bzip2 libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  && apt-get clean

# directory setup
RUN mkdir -p /opt/kpi/input /opt/kpi/raw /opt/kpi/output /opt/kpi/src /opt/bin /opt/jars

# install stuff
RUN cd /opt \
  && git clone https://github.com/droe-nmdp/kpi.git \
  && wget http://www.nextflow.io/releases/v19.01.0/nextflow?a=install -o /opt/bin/nextflow \
  && chmod 755 /opt/bin/nextflow \
  && cd /opt/bin \
  && wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz \
  && tar -zxvf KMC3.linux.tar.gz \
  && rm -f KMC3.linux.tar.gz \
  && cd /opt/jars \
  && wget http://www.apache.org/dist/commons/math/binaries/commons-math3-3.6.1-bin.tar.gz \
  && tar -zxvf commons-math3-3.6.1-bin.tar.gz \
  && rm -f /opt/jars/commons-math3-3.6.1-bin.tar.gz

# google guava
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN cd /opt \
  && git clone https://github.com/google/guava.git \
  && cd guava/guava \
  && mvn install

# env vars
ENV NXF_OPTS "-Xms4G -Xmx20G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH

# kpi source
ENV PATH /opt/bin:/opt/kpi:/opt/kpi/src:/opt/bbmap:/opt/canu/Linux-amd64/bin:$PATH
ENV CLASSPATH /opt/guava/guava/target/guava-HEAD-jre-SNAPSHOT.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:$CLASSPATH

CMD ["/opt/kpi/main.nf"]
