c     strobel.f:      
c--------------- Schumann-Runge continuum
      data fscr/1.1/,abscr/1.e-17/,escr /0.41/,
     *     eilm/0.98/,eism/0.43/,
     $     abl/2.9e-19/,abm/1.7e-18/,abk /1.15e-17/
c---------------- parameters to calculate Chapman function
     *     aa /1.0606963/,bb/0.55643831/,dd/1.7245609/,
     *     ccc/1.0619896/,f /0.56498823/,gs /0.06651874/
      DATA F_h2o/5*12.1587,6.5070,4*10.73,5*23.8226,19.2689,8*43.7116,
     $           4*35.7886,8*135.0955,6*239.2806,7*222.9263,
     $           138.789,182.3105,101.2186,72.2298/
      DATA a_h2o/0.09312,0.16262,0.14331,0.30720,0.02,1.,
     $           0.07874,0.24765,0.37943,0.27553, 
     $           0.19212,0.27513,0.24808,0.22628,0.04,1.,
     $   0.03838,0.08791,0.16378,0.11897,0.20201,0.15096,0.16933,0.004,
     $           0.01606,0.08166,0.06004,0.39433,
     $   0.04476,0.11198,0.11804,0.06353,0.08757,0.05028,0.16530,0.010,
     $           0.00532,0.02794,0.04138,0.18033,0.19552,0.18322,
     $   0.00361,0.01292,0.04997,0.08942,0.10164,0.14075,0.23042,
     $           4*1./ 
      DATA b_h2o/4.00000,0.25461,0.01489,0.00052,1000.,3.524,
     $           40.0000,2.92354,0.38207,0.04629,
     $           50.0000,4.75654,0.37403,0.02046,1000.,0.01443,
     $   100.000,6.13107,0.83251,0.23363,0.05894,0.00447,0.00036,4000.,
     $           4.00000,0.28556,0.05189,0.00629,
     $   100.000,11.2462,2.06191,0.64058,0.23003,0.08793,0.01722,500.,
     $           50.0000,7.35556,2.05710,0.40888,0.06036,0.01438,
     $     3.000,0.70579,0.15666,0.03666,0.01463,0.00534,0.00135,
     $           0.00193,0.00203,0.00004,0.00006/
      DATA p0_h2o/4*3.00,2*1000.,4*500.,4*50.,2*1000.,
     $            7*40.0,1000.,4*700.,7*1013.,1000.,6*1013.,7*700.,
     $            4*1000./ 
      DATA cm_h2o/4*0.84,2*0.00,4*0.960,4*0.6,2*0.000,
     $            7*0.88,0.000,4*0.91,7*0.800,0.000,6*0.490,7*0.38,
     $            4*0.000/ 
c     ----------------------------------------------------------------
C      data f0sun /1365./
c     ----------------------------------------------------------------

      f0sun=SOLC
      print*,'sol. const.: ',f0sun

csl   data f0sun /1365./ eventuell erhöhen für solar_max (1365 + 2 W/m²)

csl       data f0sun /1367./

csl  ---------------------------------------------------------------------

csl   2000 W /m²  extrem_test

csl      data f0sun /2000./

csl  ---------------------------------------------------------------------