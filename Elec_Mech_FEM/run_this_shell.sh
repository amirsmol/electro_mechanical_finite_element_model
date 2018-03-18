#!/bin/bash
foldername=$(date +%S_%M_%H_%d_%m_%Y)_afc_frequency_study
mkdir ../$foldername
cp -r * ../$foldername

cd ../$foldername
nohup time make all
nohup shopt -s extglob
nohup rm -rf -- !(gnuplot)
