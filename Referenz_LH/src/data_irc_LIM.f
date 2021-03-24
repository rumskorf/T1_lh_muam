c     data_irc_LIM to include in ircool_60.f

c      common /pro/ epso3(33),epsco2(70),epsh2o(kh2o)
c     *		  ,epscox(61),epscoxo(14),coolno(43)

      common /pro/ epso3(33),epsco2(70),epsh2o(25)
     *		  ,epscox(61),epscoxo(14),coolno(43)

c      common /ncirpar/indexx(11),atomo(23),dkute(23),
c     *      acox(39,12),bcox(39,12),ccox(39,12),
c     *      mina2(33),maxa2(33),minax(39),
c     *      index2(33,12),xlnp(72),
c     *      ao3(33,12),bo3(33,12)

      dimension 
     *          T118(118),T(72),Th2o(72),bl(12),cpara(12),
     *          cool(118),cooli(118),dnue(17),xnue(17)
  
      dimension tequib(kgit+1),x(kgit+1),y(kgit+1)
c
      dimension bnue(17,57),at(17),bt(17),xnc(4)

      data at /7210.3,6024.8,1614.1,139.03,21.64,2.919,
     *         0.386,0.0715,
     *         12.65,134.4,632.9,331.2,434.1,136.,35.65,9.015,1.529/
      data bt /0.182,0.094,0.081,0.08,0.068,0.06,0.059,0.067,
     *          .089,.23,.32,.296,.452,.359,.165,.104,.116/

      data coaat/718.7/
      data cobt/0.448/

c     DRUCKVERBREITERUNG 0.5-1.0 laborwerte(Goody'64)
c      data xnc /0.5,0.5,0.5,0.5/  
      data xnc /1.,1.,1.,1./  
 
 
      data dnue/120.,120.,100.,120.,100.,120.,80.,100.,
     $          150,7*100.,150./
      data xnue/100.,220.,330.,440.,550.,660.,760.,850.,
     *          1275.,1400.,1500.,1600.,1700.,1800.,1900.,2000.,2125./
      data codnue/170./
      data coxnue/667./
c      end
c
