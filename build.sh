#!/usr/bin/bash
# 
# Build file for JayCI
#
#> ./build.sh [target] 
if [ $1 == "test" ]; then
      echo ">>>>>>>>>>>> Building test.x. >>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      
      cd ./source
      make test
      cd ../
elif [ $1 == "driver1" ]; then
      echo ">>>>>>>>>>>> Building jayCI1.x. >>>>>>>>>>>>>>>>>>>>>>>>>>"
      cd ./source
      make driver1
      cd ../
elif [ $1 == "driver2" ]; then
      echo ">>>>>>>>>>>> Building jayCI2.x. >>>>>>>>>>>>>>>>>>>>>>>>>>"
      cd ./source
      make driver2
	make clean2
      cd ../
elif [ $1 == "all" ]; then
      echo ">>>>>>>>>>>> Building all: jayCI1.x, jayCI2.x >>>>>>>>>>>>"
      cd ./source
      make driver1
      make driver2
      make clean2
      cd ../
fi
echo "Build complete."
echo "Executatbles are located in ./bin diretory."
