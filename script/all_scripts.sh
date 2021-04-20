#!/bin/bash

build_folder="../build/"

ln -s /global_data/maier/common_crawl_17-04/CC-MAIN-20160205193905-00000-ip-10-236-182-209.ec2.internal.warc input

cp ${build_folder}/time/*  bin/
cp ${build_folder}/del/*   bin/
cp ${build_folder}/eps/*   bin/
cp ${build_folder}/crawl/* bin/
cp ${build_folder}/displ/* bin/
cp ${build_folder}/sing/*  bin/


./cache_test.sh
./growing.sh
./non_growing.sh
#./aggregation.sh
./deletion.sh
./load_bound.sh
