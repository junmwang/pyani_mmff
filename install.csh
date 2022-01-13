#!/bin/csh
# The source code can be compiled by using gcc 
gcc ani2ani.c -lm -o ani2ani 
gcc ani2gcrt.c -lm -o  ani2gcrt
gcc gcrt2ani.c -lm -o gcrt2ani
