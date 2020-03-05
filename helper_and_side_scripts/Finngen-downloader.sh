#!/bin/bash

# This script is meant to automatically download files from FinnGen

filestring=$(cat  FinnGen_first_download.txt |tr "\n" " ")
array=($filestring)

for i in ${array[@]};
do
wget $i
done





