#!/usr/bin/env bash

for i in Figures_ESI/S* ; do
	ln -s Code $i/Code
	rm $i/commun
	echo "$i"
done
ln -s Code Figures_Article/Code
