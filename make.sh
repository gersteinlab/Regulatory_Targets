#!/bin/bash
# -DVERSION=1.03
#make clean && make VER=1.03

g++ -g -std=c++0x -Ialglib  -Itabix-0.2.6 -Ltabix-0.2.6 -o DRM2Target DRM2Target.cpp tabixSearch.cpp gff.cpp vmstat.cpp -ltabix -lm -lz -lgsl -lgslcblas libboost_regex.a alglib.a
