FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>
# add(todo): https://github.com/google/guava
# todo: remove the supercsv stuff

# install stuff
RUN apt-get update && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip maven subversion gzip openjdk-8-jdk groovy wget \
     zlib1g-dev \
  && apt-get clean \
  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
  && cd /opt && wget https://sourceforge.net/projects/bbmap/files/latest/download \
  && mv download BBMap.tar && tar -xvzf BBMap.tar && rm BBMap.tar \
  && mkdir -p /opt/jars && cd /opt/jars \
  && wget https://github.com/google/guava/releases/download/v21.0/guava-21.0.jar \
  && apt-get clean

# env vars
ENV NXF_OPTS "-Xms1G -Xmx5G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH

# kpi source
ADD *.nf /opt/kpi/
ADD input /opt/kpi/input/
ADD output /opt/kpi/output/
ADD src /opt/kpi/src/
ENV PATH /opt/bin:$PATH
ENV PATH /opt/kpi:$PATH
ENV PATH /opt/kpi/src:$PATH
ENV PATH /opt/bbmap:$PATH
ENV CLASSPATH /opt/jars/guava-21.0.jar:$CLASSPATH

CMD ["/bin/bash"]
#CMD ["/opt/kpi/kpil.nf"]
