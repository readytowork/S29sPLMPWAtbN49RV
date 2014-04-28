#!/bin/bash
a=0
substring="Saving results to result"
fileprint=""
space=" "
pass=" 1"
notpass=" 0"
nline=$'\n'
ispassp=0
while read line
do a=$(($a+1));
  if ! [ "${line/$substring}" = "$line" ]
  then
    if [ "$ispassp" == "0" ]; then
      fileprint=$fileprint$notpass$nline
    fi
    continue
  fi
  for word in $line
  do
    if [ "$word" == "=" ]; then
      continue
    fi
    if [ "$word" == "h1" ]; then
      continue
    fi
    if [ "$word" == "h2" ]; then
      continue
    fi
    if [ "$word" == "PASS" ]; then
      fileprint=$fileprint$space$pass$space$nline
      ispassp=1
    else
      fileprint=$fileprint$word$space
      ispassp=0
    fi
  done
done < "QpskQpsk.txt"
echo "$fileprint" >> "output2.txt"
echo "Final line count is: $a";
#Q16Q16.txt qamqpsk.txt QpskQpsk.txt
