#!/bin/bash

 cat $1 |while read line
 do
   arr=($line)
   filename=${arr[0]}
   submitterid=${arr[1]}
   gunzip -c ./${filename} > ./file_upper/${submitterid}.count
 done
