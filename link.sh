#!/usr/bin/env bash

for i in Figures_ESI/S* ; do
	ln -s Code $i/Code
	echo "$i"
done
ln -s Code Figures_Article/Code
