#!/bin/bash

name=$USER
bc_year='1983' # set year for boundary conditions
year='2005' #set year for ozone and co2 conditions

datum=$(date +"20%y%m%d")
start=$(date +"20%y-%m-%d %T")
startsec=$(date +%s) #get time in sec since 1970-01-01 00:00:00 UTC
rootdir="/gemeinsam_tmp/VACILT/"
zusatz="_LH_version_10km_assimilation"
what="EL" # La ; reprecents ENSO bc for LH parameterization
outdir=${rootdir}"muam_experimental_${what}_"${name}_${bc_year}
workdir=${outdir}/${datum}${zusatz}

echo 'start: '$start
echo $datum

if [ ! -d ${workdir} ] #if folder does not already exist
then 
  if [ ! -d ${outdir} ] #if folder for me does not already exist
  then
      mkdir ${outdir} #create my own folder
  fi  
  mkdir ${workdir} #creates your own folder 
  mkdir ${workdir}/nc #subfolder for nc-data
  mkdir ${workdir}/f2 #subfolder for f2-data
  mkdir ${workdir}/dat #subfolder for other data
fi
chmod -R 750 ${workdir} #set rights for user (rwx) and group (r) including all subfolders


#./imake48.sh #execute imake to compile muam

cd ../run/ #change directory from here (mak) to run



for im in Jan #Jan #Jul #Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec #do for all selected months
do

	#define day of year for first day of each month:
	if [ $im == "Jan" ]
	then
		j=1
		mon='01'
	elif [ $im == "Feb" ]
	then
		j=32
		mon='02'
	elif [ $im == "Mar" ]
	then
		j=61
		mon='03'
	elif [ $im == "Apr" ]
	then
		j=92
		mon='04'
	elif [ $im == "May" ]
	then
		j=122
		mon='05'
	elif [ $im == "Jun" ]
	then
		j=153
		mon='06'
	elif [ $im == "Jul" ]
	then
		j=183
		mon='07'
	elif [ $im == "Aug" ]
	then
		j=214
		mon='08'
	elif [ $im == "Sep" ] 
	then 
		j=245
		mon='09'
	elif [ $im == "Oct" ]
	then
		j=275
		mon='10'
	elif [ $im == "Nov" ]
	then
		j=306
		mon='11'
	elif [ $im == "Dec" ]
	then
		j=336
		mon='12'
	fi
	
	#rm *.dx #remove old .dx-files
	
	#lower boundary conditions into run-folder - choose an option:
	#cp /projekt2/hochatm/ncep/muamfertig_alt/ncep_lb_2000-2011_$im/*.dat . #use decadal mean of NCEP data (old version PHoffmann)
	#cp /projekt2/hochatm/ncep/low_ncep_decade/ncep_2000-2010_$mon/*.dat . #use decadal mean of NCEP data (new version FLilienthal)
	#cp /projekt2/hochatm/ncep/ncep_lb_2002_$im/*.dat . #use NCEP data of a single year (2002)
	#cp /projekt2/hochatm/eraint/ERA_ENSO/EL/low_era_enso/era_1983-2010_$mon/*.dat . #use decadal mean of ERA data
	cp ${rootdir}bc_muam/${bc_year}_${mon}/*.dat . # copy paricular year bc

	cp ${rootdir}ENSO_data/*${what}*.dx . # copy ENSO bc
	for f in *.dx; do mv "${f}" "${f/_${what}/}"; done # rename ENSO bc

	rm *f2_0_*_new* #remove old f2-file
	#cp /projekt1/hochatm/MUAM_ini/"f2_0_"$im"_new"./"f2_000_"$im"_new" #copy file from ini to here
	cp ${rootdir}MUAM_ini/56-level/"f2_0_"$im"_56" ./"f2_000_"$im"_new" #use 56 level version
	
	#uncomment following line only if not starting at day 0. insert starting day
	#let oldsec=330 #270

# do for each selected section (number defines day to stop)
  for sec in 120 180 240 270 300 330 360
  do
	
	#define initializing file:
	if [ $sec == "120" ]
	then
		fold="f2_000_"$im"_new" 
	else
		fold="f2_"$oldsec"_"$im"_tst"
	fi

	#define some strings for file names
	u=_
	add1=_tst
	add2=_org	
	
	let ndek=j+14 #middle of the month
	#define new parameters for control file muam_mod_clear.txt (save in temporary files):

	#s/<text1>/<text2/ = substitute text1 by text2 ; g = global ; file1 > file2 = from source file1 to destination file2 (file1#file2)
	#-i -> in same file, -e = script commands first (not always necessary)
	cp muam_mod_clear.txt muam_mod$sec.txt #new file to change parameters
	let nend=sec*24*16 #integer calculation
	echo 'YEAR '$year
	sed -i -e 's/ NEND   = xxxxx/ NEND   = '$nend'/g' muam_mod$sec.txt
	sed -i -e 's/ NDEK   = xxx/ NDEK   = '$ndek'/g' muam_mod$sec.txt
	sed -i -e 's/ YEAR   = xxxx/ YEAR   = '$year'/g' muam_mod$sec.txt
	sed -i -e 's/ GWAM   = xxx/ GWAM   = '$im'/g' muam_mod$sec.txt
	sed -i -e 's/ FOLD   = "xxxxxxxxxxxxxx"/ FOLD   = "'$fold'"/g' muam_mod$sec.txt
	sed -i -e 's/ FNEW   = "xxxxxxxxxxxxxx"/ FNEW   = "'"f2_"$sec$u$im$add1'"/g' muam_mod$sec.txt
	sed -i -e 's/ FO3    = "o3data_xxxx.txt"/ FO3    = "'"o3data_"$year'.txt"/g' muam_mod$sec.txt
	sed -i -e 's/ FUVT   = "xxxxxxxxxxxxxxxxx"/ FUVT   = "'"uvt"$sec$u$im$add1".dx"'"/g' muam_mod$sec.txt
	if [ $sec == '300' -o $sec == '330' -o $sec == '360' ] 
	then
		sed -i -e 's/ NOUT   = xx/ NOUT   = 16/g' muam_mod$sec.txt #2h output for real simulation
	else
		sed -i -e 's/ NOUT   = xx/ NOUT   = 64/g' muam_mod$sec.txt #4h output for tuning
	fi	
			   	
	cp muam_mod$sec.txt muam_mod.txt #copy muam_mod_sec into muam_mod, ready for executing muam

	./muam48_Jan_2005 #execute muam	

	#rename retrieved files for muam2netcdf.py
	mv uvt$sec$u$im$add1.dx uvt_$im$sec.dx
	mv gwau.dx gwacu_$im$sec.dx
	mv gwav.dx gwacv_$im$sec.dx
	mv gwat.dx gwt_$im$sec.dx
	mv gwfv.dx gwfluxv_$im$sec.dx
	mv gwfu.dx gwfluxu_$im$sec.dx
	mv phi.dx  phi_$im$sec.dx
	mv wvel.dx wvel_$im$sec.dx	
	
	
	if [ $sec == '360' ]
	then
		#change pythonscript
		cp muam2netcdf_56.py muam2netcdf_temp$sec.py #create tempoary working file
  		cp muam2netcdf_56_gw.py muam2netcdf_gw_temp$sec.py #create tempoary working file
		
		sed -i -e 68c'runs = '\'$im$sec\''' muam2netcdf_temp$sec.py
		sed -i -e 68c'runs = '\'$im$sec\''' muam2netcdf_gw_temp$sec.py

		#run python for last block to create .nc files from (binary) .dx files
		python muam2netcdf_temp$sec.py
    		python muam2netcdf_gw_temp$sec.py
	fi

	#cleaning up the folder
	rm *temp*.py #remove temporary python scripts
	cp f2_$sec$u$im$add1 ${workdir}/f2/f2_$sec$u$im$add1 #move new f2-file into f2
	mv *.nc ${workdir}/nc/ #
	
	let oldsec=sec #remember previous time step

  done  
 
done

#cleaning up the folder
 rm f2_* #remove f2-files
 rm *.dat #remove old .dat files in folder (lower boundary conditions)
 rm *.dx

ende=$(date +"20%y-%m-%d %T") #get Date and Time
endsec=$(date +%s) #get time in sec since 1970-01-01 00:00:00 UTC
echo "stop: "$ende", "$endsec

#calculate time for the run:
time=$(( $endsec - $startsec ))
let timeh=time/3600 #in hours
let timem=(time-timeh*3600)/60 #in minutes
let times=(time-timeh*3600-timem*60) #in seconds
echo $timeh" "$timem" "$times

