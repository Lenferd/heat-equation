#!/bin/bash 
PATH_PROGRAM=$(find ../../build -name 'app_Sergey_OpenMP_Euler')
echo "$PATH_PROGRAM"
cp -f "$PATH_PROGRAM" .

for (( h = 8; h <= 1024; h *= 2 ))
do
	echo "*********"
	echo $h
	echo "*********"
	sed -i 's/NX=.*/NX='"$h"'/g' setting2.ini
	sed -i 's/NY=.*/NY='"$h"'/g' setting2.ini
	sed -i 's/NZ=.*/NZ='"$h"'/g' setting2.ini
	python3 generation2.py

	for (( dt = 3; dt <= 4; dt++ ))
	do
		echo "-------"
		echo $dt
		echo "-------"
		sed -i 's/dt=.*/dt=1e-'"$dt"'/g' setting2.ini
		./app_Sergey_OpenMP_Euler
	done
done
rm res.txt
