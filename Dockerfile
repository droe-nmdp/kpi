FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>
# todo: change the supercsv stuff

# install stuff
# https://github.com/gmarcais/Jellyfish
RUN apt-get update && apt-get install -qyy curl git make cmake \
     gcc g++ unzip maven subversion gzip openjdk-8-jdk groovy wget \
  && apt-get clean \
  && cd /opt && mvn dependency:get -Dartifact=net.sf.supercsv:super-csv:2.4.0 \
  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
  && cd /opt/ && wget http://www.cs.cmu.edu/~ckingsf/software/bloomtree/sbt-binary.tar.gz \
  && tar -zxvf sbt-binary.tar.gz \
  && ln -s /opt/sbt-binary/bt /opt/bin \
  && ln -s /opt/sbt-binary/get_bfsize.sh /opt/bin \
  && rm /opt/sbt-binary.tar.gz \
  && apt-get clean

RUN git clone https://github.com/droe-nmdp/kpi /opt/kpi/repos

RUN ln -s /opt/kpi/repos/src /opt/kpi/
RUN cp /opt/kpi/repos/*.nf /opt/kpi
RUN rsync -a /opt/kpi/repos/input /opt/kpi/

# env vars
ENV NXF_OPTS "-Xms1G -Xmx5G"

# kpi source
ENV PATH /opt/bin:$PATH
ENV PATH /opt/kpi:$PATH
ENV PATH /opt/kpi/src:$PATH
ENV CLASSPATH /opt/kpi/src/super-csv.jar:$CLASSPATH

CMD ["/bin/bash"]
#CMD ["/opt/kpi/kpil.nf"]
