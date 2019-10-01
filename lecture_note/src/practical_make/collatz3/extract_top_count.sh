#!/bin/sh

inputfile=$1

cat $inputfile | awk '{print $2, $3}' | sort -nr | head -1
