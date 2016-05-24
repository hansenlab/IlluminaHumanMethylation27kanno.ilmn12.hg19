#!/bin/bash -e
cd extdata

curl -O ftp://hgdownload.soe.ucsc.edu/apache/htdocs/goldenPath/hg19/database/snp146Common.txt.gz
curl -O ftp://hgdownload.soe.ucsc.edu/apache/htdocs/goldenPath/hg19/database/snp144Common.txt.gz
curl -O ftp://hgdownload.soe.ucsc.edu/apache/htdocs/goldenPath/hg19/database/snp142Common.txt.gz
curl -O ftp://hgdownload.soe.ucsc.edu/apache/htdocs/goldenPath/hg19/database/snp141Common.txt.gz
curl -O ftp://hgdownload.soe.ucsc.edu/apache/htdocs/goldenPath/hg19/database/snp138Common.txt.gz

gunzip -c snp146Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp146Common_small.txt.gz
gunzip -c snp144Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp144Common_small.txt.gz
gunzip -c snp142Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp142Common_small.txt.gz
gunzip -c snp141Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp141Common_small.txt.gz
gunzip -c snp138Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp138Common_small.txt.gz






