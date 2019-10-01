#!/bin/sh

inputfile=$1

cat $inputfile | sort -nr | head -1
