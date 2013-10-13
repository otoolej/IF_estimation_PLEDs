#!/bin/bash

LINEN=`grep -n '^function' $1 | head -n 1 | cut -f1 -d:`
echo $LINEN

sed '\${LINEN}i license_header.m' $1

