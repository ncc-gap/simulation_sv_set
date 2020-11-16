FROM python:3.7.9-buster

MAINTAINER ken-ichi chiba <kchiba@hgc.jp>

RUN apt-get -y update
RUN apt-get install -y wget

WORKDIR /usr/local

RUN wget https://github.com/fritzsedlazeck/SURVIVOR/archive/v1.0.6.tar.gz && \
    tar xzvf v1.0.6.tar.gz && \
    cd SURVIVOR-1.0.6/Debug && \
    make

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && \
    tar jxvf minimap2-2.17.tar.bz2 && \
    cd minimap2-2.17 && \
    make

RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar jxvf samtools-1.11.tar.bz2 && \
    cd samtools-1.11 && \
    ./configure && make && make install

RUN wget https://github.com/bcgsc/NanoSim/archive/v2.6.0.tar.gz && \
    tar xzvf v2.6.0.tar.gz && \
    cd NanoSim-2.6.0/src && \
    grep -l 'from sklearn.externals import joblib' *py  | xargs sed -i.bak -e 's;from sklearn.externals import joblib;import joblib;g'

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN pip install six==1.15.0
RUN pip install pysam==0.15.4
RUN pip install numpy==1.18.5
RUN pip install scipy==1.4.1
RUN pip install HTSeq==0.10.0
RUN pip install scikit-learn==0.23.2

ENV PATH $PATH:/usr/local/SURVIVOR-1.0.6/Debug
ENV PATH $PATH:/usr/local/minimap2-2.17
ENV PATH $PATH:/usr/local/samtools-1.9/bin


