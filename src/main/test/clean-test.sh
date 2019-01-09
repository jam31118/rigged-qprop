rm -rf ./ati-ppp
make clean-data -f ../../makefile -C ./ati
make clean-data -f ../../makefile -C ./ati-before-ppp
make clean -C ..
rm ../ppp
