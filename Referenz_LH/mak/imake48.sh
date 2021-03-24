cd ../src/

'rm' *.o
'rm' muamlib.a

export PATH=${PATH}:/opt/intel/Compiler/11.1/056/bin/intel64/
ifort -fpe0 -c -g -traceback -O2 -xHost acc_ASM.f bi_diff_60.f CCO2GR.f egwd_RSHU.f fourier_filter.f fourier_filter_lb.f fourier_matrix.f geopot_60.f grwaves_RSHU.f heating_LH.f init_60.f irc_60.f ircool_60.f legandr_philb.f levy.f model_48.f molcon_60_1.f ncool_PWs_60.f O3model_Berlin_60.f open_60.f pcoolg.f plwaves_LIM.f species_60.f strobel_60.f tendenz_60_1.f tendenz_60_a.f vbew.f 
#ifort -fpe0 -c bi_diff_60.f CCO2GR.f fourier_filter.f fourier_filter_lb.f fourier_matrix.f geopot_60.f grwaves_RSHU.f init_60.f irc_60.f ircool_60.f legandr_philb.f levy.f model_48.f molcon_60_1.f ncool_PWs_60.f O3model_Berlin_60.f open_60.f pcoolg.f plwaves_LIM.f species_60.f strobel_60.f tendenz_60_1.f tendenz_60_a.f vbew.f 

ar -q muamlib.a

rm model_48.o

ar -r muamlib.a *.o

ifort -fpe0 -g -traceback -O2 -xHost -o muam48_Jan_2005 model_48.f muamlib.a sublib.a
#ifort -fpe0 -o muam48 model_48.f muamlib.a sublib.a
 
mv muam48_Jan_2005 ../run/

#export PATH="${PATH}:/home/hoffmann/tools/bin:/opt/intel/Compiler/11.1/056/bin/intel64"
