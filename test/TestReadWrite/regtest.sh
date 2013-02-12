#!/bin/bash

myApp=$1
testDir=`dirname $2`

cd /tmp
readarray < $2 &> LR_regtest.log
$myApp $MAPFILE 

diff -u TestReadWrite.lr TestReadWrite2.lr
result=$?
if [ $result -eq 0 ]; then
	diff -u TestReadWrite.lr TestReadWrite3.lr
	result=$?
fi

rm -f TestReadWrite.lr TestReadWrite2.lr TestReadWrite3.lr TestReadWrite4.lr

test $result -eq 0 && exit 0
exit 1
