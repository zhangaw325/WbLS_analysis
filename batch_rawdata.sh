#!/bin/bash

count=0
for name in `cat ./filenamelist_runs20190408.dat`
do
    #echo $name
    thisarray[$count]=$name
    ((count = count + 1))
done
   /gpfs01/lbne/users/aw325/sw/Python2.7.13.install/bin/python2.7 Read_MultipleFiles_WriteTree_v9.py ${thisarray[${1}]}
