FROM python:3.7.9-buster

MAINTAINER ken-ichi chiba <kchiba@hgc.jp>

RUN apt-get -y update
RUN apt-get install -y wget

WORKDIR /usr/local

RUN mkdir -p ucsc_tools && \
    cd ucsc_tools && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit && \
    chmod +x faToTwoBit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo && \
    chmod +x twoBitInfo && \
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod +x blat && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit && \
    chmod +x faSplit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftUp && \
    chmod +x liftUp && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/axtChain && \
    chmod +x axtChain && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainMergeSort && \
    chmod +x chainMergeSort && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainSplit && \
    chmod +x chainSplit && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/chainNet && \
    chmod +x chainNet && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/netChainSubset && \
    chmod +x netChainSubset && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod +x liftOver

ENV PATH $PATH:/usr/local/ucsc_tools

