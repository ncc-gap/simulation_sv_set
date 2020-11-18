FROM python:3.7.9-buster

MAINTAINER ken-ichi chiba <kchiba@hgc.jp>

RUN apt-get -y update
RUN apt-get install -y wget

WORKDIR /usr/local

RUN mkdir -p ucsc_tools && \
    cd ucsc_tools && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit && \
    chmod u+x faToTwoBit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo && \
    chmod u+x twoBitInfo && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod u+x blat && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit && \
    chmod u+x faSplit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftUp && \
    chmod u+x liftUp && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtChain && \
    chmod u+x axtChain && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort && \
    chmod u+x chainMergeSort && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSplit && \
    chmod u+x chainSplit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainNet && \
    chmod u+x chainNet && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netChainSubset && \
    chmod u+x netChainSubset && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod u+x liftOver

ENV PATH $PATH:/usr/local/ucsc_tools

