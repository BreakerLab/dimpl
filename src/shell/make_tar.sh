#/bin/bash

name=`pwd |sed 's#.*/##'`
cd ..
tar czf ${name}.done.tar.gz ${name}
