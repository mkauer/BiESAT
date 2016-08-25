#!/bin/sh

parnames="exp0 exp1 Cn500 Ce500 Cs500 Css500 CGn500 CGe500 CGs500 Cn1000 Ce1000 Cs1000 Css1000 CGn1000 CGe1000 CGs1000 Gn482 Ge482 Gs482 Gc482 Gn976 Ge976 Gs976 Gc976 offset"

file=par_changes.dat

i=0

for parname in $parnames
do
    #echo -e "par[$i]  -->  $parname "
    sed -i -c "s/par\[$i\]/$parname/" $file
    #sed -i -c "s/par\[$i\]/$parname/" $1
    ((i++))
done


