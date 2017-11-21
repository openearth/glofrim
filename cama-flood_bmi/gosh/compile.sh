#!/bin/sh 

clean=$1       # 1th argument in script : if == "yes" then clean before compiling 

#CAMADIR="/Users/yamadai/work/CaMa-Flood"
CAMADIR=`pwd`/..
BASE=$CAMADIR  # code location

libs="mod lib src map out"
for lib in $libs
do 
  cd $BASE/$lib
  if [ $clean = "yes" ];then
    make clean
  fi
  echo "*********** $lib **********"
  make all
done 

cd $BASE/src/
make MAIN_day
if [ -r MAIN_day ]; then
 echo "Compilation OK!  Executable created: $BASE/src/MAIN_day"
else
 echo "Problems during Compilation. Executable not created."
fi

