#!/bin/bash

myApp=$1
testDir=`dirname $2`

folder=`mktemp -d -t LR_regtestXXXXX`
cd $folder
readarray < $2 &> log
$myApp $MAPFILE 

diff -u TestReadWrite.lr TestReadWrite2.lr
result=$?
if [ $result -eq 0 ]; then
	diff -u TestReadWrite.lr TestReadWrite3.lr
	result=$?
fi

rm -Rf $folder
test $result -eq 0 && exit 0

exit 1
