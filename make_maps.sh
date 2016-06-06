#!/bin/bash
# What this script does
# (1) creates the tablemaking and mapmaking executable from the fortran codes. These codes are in src/pks2map_src
# (2) uses the table making code to create a table for pasting on the halo profile. This table can be found in tables/
# (3) uses the mapmaking code to read in this table file and a halo catalogue and make 2 maps - one flatsky and one fullsky. flatsky maps will end with ___fs.map and fullsky will end with ____hp.fits. These maps can be found in mapsout/ 
# (4) I also included some plotting scripts in plotting/ 
for freq in `echo 201`
do
    for alpha in `echo 0`
    do
        #SET UP DIRECTORIES
	dir_tab=tables
	dir_out=mapsout
	dir_bin=bin
	dir_src=src
	if [ ! -d $dir_tab  ]; then mkdir $dir_tab ; fi
	if [ ! -d $dir_out  ]; then mkdir $dir_out ; fi
	if [ ! -d $dir_bin  ]; then mkdir $dir_bin ; fi
	
        #MAKE EXECUTABLE FILES
	if [ ! -f ${dir_bin}/pks2map  ]; then cd $dir_src; make pks2map ; cd ../ ; fi
	if [ ! -f ${dir_bin}/make_maptable  ]; then cd $dir_src; make make_maptable ; cd ../ ; fi
    
        #MAKE TABLE FILE FOR SPECIFIED FREQUENCY AND ALPHA
	tabfile=${dir_tab}/tab_${freq}GHz_alpha_${alpha}
	if [ ! -f $tabfile  ]; then
	    ./bin/make_maptable $tabfile 1 $freq $alpha 
	fi
    
        #MAKE MAP 
	#USAGE ./bin/pks2map <filein> <fileout> <tablefile> 
	#      [<zmin> <nside> <scramble> <center> <npix> <fov> 
	#      <zmax> <chihview> <model> <profile>]
	
	ppfile=8Gpc_n4096_nb23_nt18_merge.pksc.13579
	outfile=${dir_out}/13579_${freq}GHz_alpha_${alpha}
	if [ ! -f ${outfile}_hp.fits  ]; then
	    mpirun -np 8 ./bin/pks2map $ppfile $outfile $tabfile 0.0 4096 0 1 4096 10 1.245 0 1 1
	fi
    done
done


