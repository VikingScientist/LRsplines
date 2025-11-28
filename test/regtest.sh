#!/bin/bash

mysim=$1
BINARY_DIR=$3

cd `dirname $2`
readarray < $2
templog=`mktemp -t lrsplinelogXXXXX`
$mysim $MAPFILE 2>&1 > $templog
globres=1
IFS=$'\n'
i=0
for line in `cat $2`
do
  test -z "$line" && continue
  ((i++))
  test $i = 1 && continue # skip first line of regression file
  result=0
  if grep -q "$line" $templog
    then result=1 
  fi
  if test $result -eq 0 
  then 
    if test $globres -eq 1
    then
      echo "-------- $2 --------" >> $BINARY_DIR/failed.log
    fi
    globres=0
    echo Failed to find output: $line >> $BINARY_DIR/failed.log
  fi
done

if test $globres -eq 0
then
  cat $templog >> $BINARY_DIR/failed.log
  rm -f $templog
  exit 1
fi

rm -f $templog
exit 0
