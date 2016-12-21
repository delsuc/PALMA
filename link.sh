#!/usr/bin/env bash

for i in Figures_ESI/S* ; do
	( cd $i; ln -s ../../Code Code )
	echo "$i"
done

( cd Figures_Article; ln -s ../Code Code)
