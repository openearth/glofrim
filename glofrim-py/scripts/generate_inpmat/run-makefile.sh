#!/bin/sh
MKINCLUDE=$1
LF=$(printf '\\\012_')
LF=${LF%_}
rm -f Makefile
sed -e '2d' sample_Makefile > Makefile
sed -i -e "2s%^%include   ${MKINCLUDE}${LF}%" Makefile
rm -f Makefile-e
make clean
make all
