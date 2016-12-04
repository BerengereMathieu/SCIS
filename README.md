# SCIS
Author Bérengère MATHIEU

## Licence 
This program is free software; you can redistribute it and/or modify
it under the terms of the CeCILL-C FREE SOFTWARE LICENSE  as published by
here http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
CeCILL-C FREE SOFTWARE LICENSE  for more details.

You should have received a copy of the CeCILL-C FREE SOFTWARE LICENSE
along with this program.


## Dependencies : 

	* gimp development library


### Linux and Mac

	Use gimp and gimptool to automatically compile and install SCIS :

	* for linux users : install package gimp  and libgimp2.0-dev 
	* for mac users : you need to install gimp with a packages manadger, like Macport
	* CC=gcc CFLAGS="-v -ansi" LIBS=-lm gimptool-2.0 --install SCIS.c

### Windows
	
	* install gimp software 
	* copy SCIS_gimp_plugin_win32 file in the gimp plugin folder. Gimp main folder is usually somewhere in Program Files. Once in the GIMP's main folder navigate to lib\gimp\*version*\ where as *version* represents the version of The Gimp. You should see a directory named **plug-in**. You can also try to locate a well known Gimp plugin like **file-jpeg** plugin. 

## Use
	* open the image with gimp
	* create a new layer 
	* give seeds to algorithm by drawing strokes on this layer : 
		* use one color by semantic class ;
		* be careful to use uniform color and avoid drawing tools with transparency effect (you can use paintbrush or pencil with brush **2. Hardness 100**)
	* run SCIS : filter -> segmentation -> SCIS
	
	

## Acknowledgement 


I thank Pedro Felzenszwalb of Brown University, for allowing me to use its
implementation of its over-segmentation algorithm (http://cs.brown.edu/~pff/segment/)

I also gratefully acknowledge Chih-Chung Chang and Chih-Jen Lin of
National Taiwan University, for allowing me to use their SVM training
and classification library : libSVM (https://www.csie.ntu.edu.tw/~cjlin/libsvm/).
