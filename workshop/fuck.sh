#!/bin/bash
for ((i=5;i<=15;i++));
do
   for ((j=i;j<=15;j++));
   do
      echo $((i*100))" "$((j*100))
      ./exp06 524288 $((i*100)) $((j*100))
   done
done
