# Ranking neuroblastoma copy-number altered genes

## Get started

    docker build -t nbranking .
    docker run -v /c/Users/christophe/data/FilesChristophe:/root/.local/share/lostdata/private -it nbranking /bin/bash

## For development

    docker run -v \
    /c/Users/christophe/data/FilesChristophe:/root/LSData/ -v \
    /c/Users/christophe/repos:/code -it nbranking /bin/bash
