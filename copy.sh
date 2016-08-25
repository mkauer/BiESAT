#!/bin/bash

if [ -z $1 ];then
    echo "need a destination dir"
    exit
fi

path="$1"
#mydate=`date +%y.%m.%d`
mydate='09.07.19'
dir="${mydate}_BiESAT"

if [ -d "$path" ];then
    mkdir -p "$path/$dir"
    cp -vf AA_* "$path/$dir"
fi

