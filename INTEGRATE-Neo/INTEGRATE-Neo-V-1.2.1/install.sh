#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo
    echo "        ./intall.sh -o destination_dir"
    echo
    echo "        Options :"
    echo "        -o: string   [the object directory of installing INTEGRATE-Neo]"
    echo
    exit
fi

args=("$@")

destination_dir=${args[1]}

echo "making destination dir..."

if [ ! -d $destination_dir/tmp ]; then
    mkdir -p $destination_dir"/tmp"
fi

echo "copying files..."

cp -r src/* $destination_dir
cd $destination_dir

echo "Compiling "

cd tmp
mv ../BedpeAnnotator ./
mkdir BedpeAnnotator_build
cd BedpeAnnotator_build
cmake ../BedpeAnnotator -DCMAKE_BUILD_TYPE=release
make

cd ..
mv ../Subsetter ./
mkdir Subsetter_build
cd Subsetter_build
cmake ../Subsetter -DCMAKE_BUILD_TYPE=release
make

cd ..
cp BedpeAnnotator_build/bin/fusionBedpeAnnotator ../
cp Subsetter_build/bin/fusionBedpeSubsetter ../

cd ..

echo "adding +x mode"
chmod +x HLAminerToTsv.py
chmod +x fusionBedpeAnnotator
chmod +x fusionBedpeSubsetter
chmod +x integrate-neo.py
chmod +x runAddNetMHC4Result.py
chmod +x runHLAminer.py
chmod +x runNetMHC4WithSMCRNABedpe.py

echo "removing temporary files"

rm -rf tmp

echo 
echo 
echo "Please edit the setup.ini file"
echo "Enjoy!"
echo
echo 


