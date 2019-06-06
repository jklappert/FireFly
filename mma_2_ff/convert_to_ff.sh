#!/bin/bash

mkdir -p ff_conv

funs=$1
vars=$2
threads=$3

math -run "functions=\"$funs\";variables=\"$vars\";nthreads=$threads" -script convert_to_ff.m

for f in ff_conv/*.cpp; do
    echo "Converting $f"
    sed -i 's/"mpz_class(",/mpz_class(/g' $f
    sed -i 's/,")"/)/g' $f
done

