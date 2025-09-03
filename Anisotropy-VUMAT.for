      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
      parameter( zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0)
      parameter(third = one/three, half = .5d0, twothirds = two/three)
      parameter(threehalfs = 1.5d0, sqrt23=0.816496580927726032732d0)
      parameter(tol=1.d-6)
C

      
C     Input elastic properties
      e1     = props(1)
      e2     = e1
      e3     = props(2)
      pr12   = props(3)
      pr13   = props(4)
      pr23   = pr13
      g12    = props(5)
      g31    = props(6)
      g23    = g31
C     Compute reciprocal Poisson＊s ratios to ensure symmetry
      pr21 = pr12 * e2 / e1
      pr31 = pr13 * e3 / e1
      pr32 = pr23 * e3 / e2
C     Compute denominator for stiffness matrix
      denom = one - pr12*pr21 - pr13*pr31 - pr23*pr32 
     1 - two*pr12*pr23*pr31
C     Build stiffness matrix components
      
      c11 = e1 * (one  - pr23*pr32) / denom
      c22 = e2 * (one  - pr13*pr31) / denom
      c33 = e3 * (one  - pr12*pr21) / denom
      c12 = e1 * (pr21 + pr31*pr23) / denom
      c13 = e1 * (pr31 + pr21*pr32) / denom
      c23 = e2 * (pr32 + pr31*pr12) / denom
      c21 = c12
      c31 = c13
      c32 = c23
C     Shear stiffness directly from Gij
      c44 = g12*two
      c55 = g23*two
      c66 = g31*two
C

C     Reference stress-strain
C     Setup stiffness and strain for bilinear curve, 0 degree
      stressy1       = props(7)
      stif11         = props(8)
      stif21         = props(9)
      stif31         = props(10)
      stif41         = props(11)
      stif51         = props(12)
      stif61         = props(13)
      stif71         = props(14)
      stif81         = props(15)
      stif91         = props(16)
      stif101        = props(17)
      stif111        = props(18)
      stif121        = props(19)
      stif131        = props(20)
      stif141        = props(21)
      stif151        = props(22)
      stif161        = props(23)
      stif171        = props(24) 
      stif181        = props(25)
      stif191        = props(26)
      stif201        = props(27)  
      stif211        = props(28)
      stif221        = props(29)
      stif231        = props(30)
      
      strain1        = 0.0002d0
      strain2        = 0.0005d0
      strain3        = 0.0008d0
      strain4        = 0.0011d0
      strain5        = 0.0015d0
      strain6        = 0.002d0
      strain7        = 0.004d0
      strain8        = 0.007d0
      strain9        = 0.008d0
      strain10       = 0.01d0
      strain11       = 0.02d0
      strain12       = 0.03d0
      strain13       = 0.04d0
      strain14       = 0.05d0
      strain15       = 0.07d0
      strain16       = 0.1d0
      strain17       = 0.150d0
      strain18       = 0.209081d0
      strain19       = 0.238228d0
      strain20       = 0.3d0
      strain21       = 0.4d0
      strain22       = 0.6d0
      strain23       = 0.8d0
C     Setup stiffness and strain for bilinear curve, 90 degree
      stressy2       = props(31)
      stif12         = props(32)
      stif22         = props(33)
      stif32         = props(34)
      stif42         = props(35)
      stif52         = props(36)
      stif62         = props(37)
      stif72         = props(38)
      stif82         = props(39)
      stif92         = props(40)
      stif102        = props(41)
      stif112        = props(42)
      stif122        = props(43)
      stif132        = props(44)
      stif142        = props(45)
      stif152        = props(46)
      stif162        = props(47)
      stif172        = props(48) 
      stif182        = props(49)
      stif192        = props(50)
      stif202        = props(51)  
      stif212        = props(52)
      stif222        = props(53)
      stif232        = props(54)
C     Setup stiffness and strain for bilinear curve, 45 degree
      stressy3       = props(55)
      stif13         = props(56)
      stif23         = props(57)
      stif33         = props(58)
      stif43         = props(59)
      stif53         = props(60)
      stif63         = props(61)
      stif73         = props(62)
      stif83         = props(63)
      stif93         = props(64)    
      stif103        = props(65)
      stif113        = props(66)
      stif123        = props(67)
      stif133        = props(68)
      stif143        = props(69)
      stif153        = props(70)
      stif163        = props(71)
      stif173        = props(72) 
      stif183        = props(73)
      stif193        = props(74)
      stif203        = props(75)  
      stif213        = props(76)
      stif223        = props(77)
      stif233        = props(78)
C -----------------------------------------------------------   
     
      
      do 100 i = 1,nblock
C     Read value from the last increment

      ep_old      = stateOld (i,1) 

C     Apply elastic predictor
C     忖考 = C : 忖汍 (Voigt notation)        
      stressInc11 = c11*strainInc(i,1) + c12*strainInc(i,2) 
     1              + c13*strainInc(i,3)
      stressInc22 = c21*strainInc(i,1) + c22*strainInc(i,2) 
     1              + c23*strainInc(i,3)
      stressInc33 = c31*strainInc(i,1) + c32*strainInc(i,2) 
     1              + c33*strainInc(i,3)
      stressInc12 = c44*strainInc(i,4)
      stressInc23 = c55*strainInc(i,5)
      stressInc31 = c66*strainInc(i,6)

C     考 = 考 + 忖考          
      stress11 = stressOld(i,1) + stressInc11
      stress22 = stressOld(i,2) + stressInc22
      stress33 = stressOld(i,3) + stressInc33
      stress12 = stressOld(i,4) + stressInc12
      stress23 = stressOld(i,5) + stressInc23
      stress31 = stressOld(i,6) + stressInc31
      
C -----------------------------------------------------------   
C     Hill 1948 yield function

C     Calculate isotropic yield stress
      if (ep_old.le.zero) then
          yield_old = stressy1
      else

          if (ep_old.le.strain1) then
               yield_old = stressy1 + stif11*ep_old
          else
              if (ep_old.le.strain2) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(ep_old-strain1)
              else
                  if (ep_old.le.strain3) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(ep_old-strain2)
                  else 
                      if (ep_old.le.strain4) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(ep_old-strain3)
                      else
                          if (ep_old.le.strain5) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(ep_old-strain4)
                          else
                              if (ep_old.le.strain6) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(ep_old-strain5)
                              else
                                  if (ep_old.le.strain7) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(ep_old-strain6)
                                  else
                                      if (ep_old.le.strain8) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(ep_old-strain7)
                                      else
                                          if (ep_old.le.strain9) then
                 yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(ep_old-strain8)
                                          else            
                 if (ep_old.le.strain10) then
                  yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(ep_old-strain9)
                 else
                     
                  if (ep_old.le.strain11) then
                  yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(ep_old-strain10)
                  else
                      if (ep_old.le.strain12) then
                  yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(ep_old-strain11)
                      else
                          if (ep_old.le.strain13) then
                  yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(ep_old-strain12)    
                          else
                              if (ep_old.le.strain14) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(ep_old-strain13)
                          else
                              if (ep_old.le.strain15) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(ep_old-strain14)
                    else
                              if (ep_old.le.strain16) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(ep_old-strain15)
                    else
                              if (ep_old.le.strain17) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(ep_old-strain16)              
                    else
                              if (ep_old.le.strain18) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(ep_old-strain17)  
                    else
                              if (ep_old.le.strain19) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(strain18-strain17)
     9                       + stif191*(ep_old-strain18)
                    
                    else
                              if (ep_old.le.strain20) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(strain18-strain17)
     9                       + stif191*(strain19-strain18)
     1                       + stif201*(ep_old-strain19)
                   
                    else
                              if (ep_old.le.strain21) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(strain18-strain17)
     9                       + stif191*(strain19-strain18)
     1                       + stif201*(strain20-strain19)
     2                       + stif211*(ep_old-strain20)               
                    else
                              if (ep_old.le.strain22) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(strain18-strain17)
     9                       + stif191*(strain19-strain18)
     1                       + stif201*(strain20-strain19)
     2                       + stif211*(strain21-strain20) 
     3                       + stif221*(ep_old-strain21) 
                    
                    else
                              if (ep_old.le.strain23) then
                    yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(strain18-strain17)
     9                       + stif191*(strain19-strain18)
     1                       + stif201*(strain20-strain19)
     2                       + stif211*(strain21-strain20) 
     3                       + stif221*(strain22-strain21)
     4                       + stif231*(ep_old-strain22)  
                              else
                      yield_old = stressy1 + stif11*strain1 
     1                       + stif21*(strain2-strain1)
     2                       + stif31*(strain3-strain2)
     3                       + stif41*(strain4-strain3)
     4                       + stif51*(strain5-strain4)
     5                       + stif61*(strain6-strain5)
     6                       + stif71*(strain7-strain6)
     7                       + stif81*(strain8-strain7)
     8                       + stif91*(strain9-strain8)
     9                       + stif101*(strain10-strain9)
     1                       + stif111*(strain11-strain10)
     2                       + stif121*(strain12-strain11)
     3                       + stif131*(strain13-strain12)
     4                       + stif141*(strain14-strain13)
     5                       + stif151*(strain15-strain14)
     6                       + stif161*(strain16-strain15)
     7                       + stif171*(strain17-strain16)
     8                       + stif181*(strain18-strain17)
     9                       + stif191*(strain19-strain18)
     1                       + stif201*(strain20-strain19)
     2                       + stif211*(strain21-strain20) 
     3                       + stif221*(strain22-strain21)
     4                       + stif231*(strain23-strain22)
     5                       + stif231*(ep_old-strain23)          
                      endif
                  endif
                    endif
                    endif
                    endif
                     endif
                      endif
                    endif
                    endif
                    endif
                    endif
                    endif
                 endif
                 
                                          endif
                                      endif
                                  endif
                              endif
                          endif
                      endif
                  endif              
              endif
          endif
          endif
          endif
C Calculate rij in Hill model
      if (ep_old.le.zero) then
          yield2 = stressy2
      else
          if (ep_old.le.strain1) then
               yield2 = stressy2 + stif12*ep_old
          else
              if (ep_old.le.strain2) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(ep_old-strain1)
              else
                  if (ep_old.le.strain3) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(ep_old-strain2)
                  else 
                      if (ep_old.le.strain4) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(ep_old-strain3)
                      else
                          if (ep_old.le.strain5) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(ep_old-strain4)
                          else
                              if (ep_old.le.strain6) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(ep_old-strain5)
                              else
                                  if (ep_old.le.strain7) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(ep_old-strain6)
                                  else
                                      if (ep_old.le.strain8) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(ep_old-strain7)
                                      else
                                          if (ep_old.le.strain9) then
                 yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(ep_old-strain8)
                                          else            
                if (ep_old.le.strain10) then
                  yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(ep_old-strain9)
                 else
                     
                  if (ep_old.le.strain11) then
                  yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(ep_old-strain10)
                  else
                      if (ep_old.le.strain12) then
                  yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(ep_old-strain11)
                      else
                          if (ep_old.le.strain13) then
                  yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(ep_old-strain12)    
                          else
                              if (ep_old.le.strain14) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(ep_old-strain13)
                          else
                              if (ep_old.le.strain15) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(ep_old-strain14)
                    else
                              if (ep_old.le.strain16) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(ep_old-strain15)
                    else
                              if (ep_old.le.strain17) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(ep_old-strain16)              
                    else
                              if (ep_old.le.strain18) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(ep_old-strain17)  
                    else
                              if (ep_old.le.strain19) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(strain18-strain17)
     9                       + stif192*(ep_old-strain18)
                    
                    else
                              if (ep_old.le.strain20) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(strain18-strain17)
     9                       + stif192*(strain19-strain18)
     1                       + stif202*(ep_old-strain19)
                  
                    else
                              if (ep_old.le.strain21) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(strain18-strain17)
     9                       + stif192*(strain19-strain18)
     1                       + stif202*(strain20-strain19)
     2                       + stif212*(ep_old-strain20)               
                    else
                              if (ep_old.le.strain22) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(strain18-strain17)
     9                       + stif192*(strain19-strain18)
     1                       + stif202*(strain20-strain19)
     2                       + stif212*(strain21-strain20) 
     3                       + stif222*(ep_old-strain21) 
                  
                    else
                              if (ep_old.le.strain23) then
                    yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(strain18-strain17)
     9                       + stif192*(strain19-strain18)
     1                       + stif202*(strain20-strain19)
     2                       + stif212*(strain21-strain20) 
     3                       + stif222*(strain22-strain21)
     4                       + stif232*(ep_old-strain22)  
                              else
                      yield2 = stressy2 + stif12*strain1 
     1                       + stif22*(strain2-strain1)
     2                       + stif32*(strain3-strain2)
     3                       + stif42*(strain4-strain3)
     4                       + stif52*(strain5-strain4)
     5                       + stif62*(strain6-strain5)
     6                       + stif72*(strain7-strain6)
     7                       + stif82*(strain8-strain7)
     8                       + stif92*(strain9-strain8)
     9                       + stif102*(strain10-strain9)
     1                       + stif112*(strain11-strain10)
     2                       + stif122*(strain12-strain11)
     3                       + stif132*(strain13-strain12)
     4                       + stif142*(strain14-strain13)
     5                       + stif152*(strain15-strain14)
     6                       + stif162*(strain16-strain15)
     7                       + stif172*(strain17-strain16)
     8                       + stif182*(strain18-strain17)
     9                       + stif192*(strain19-strain18)
     1                       + stif202*(strain20-strain19)
     2                       + stif212*(strain21-strain20) 
     3                       + stif222*(strain22-strain21)
     4                       + stif232*(strain23-strain22)
     5                       + stif232*(ep_old-strain23)           
                              endif
                  endif
                    endif
                    endif
                    endif
                     endif
                      endif
                    endif
                    endif
                    endif
                    endif
                    endif
                 endif
                 
                                          endif
                                      endif
                                  endif
                              endif
                          endif
                      endif
                  endif              
              endif
          endif
          endif
          endif
      
      if (ep_old.le.zero) then
          yield3 = stressy3
      else
          if (ep_old.le.strain1) then
               yield3 = stressy3 + stif13*ep_old
          else
              if (ep_old.le.strain2) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(ep_old-strain1)
              else
                  if (ep_old.le.strain3) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(ep_old-strain2)
                  else 
                      if (ep_old.le.strain4) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(ep_old-strain3)
                      else
                          if (ep_old.le.strain5) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(ep_old-strain4)
                          else
                              if (ep_old.le.strain6) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(ep_old-strain5)
                              else
                                  if (ep_old.le.strain7) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(ep_old-strain6)
                                  else
                                      if (ep_old.le.strain8) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(ep_old-strain7)
                                      else
                                          if (ep_old.le.strain9) then
                 yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(ep_old-strain8)
                                          else            
                 
                  if (ep_old.le.strain10) then
                  yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(ep_old-strain9)
                 else
                     
                  if (ep_old.le.strain11) then
                  yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(ep_old-strain10)
                  else
                      if (ep_old.le.strain12) then
                  yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(ep_old-strain11)
                      else
                          if (ep_old.le.strain13) then
                  yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(ep_old-strain12)    
                          else
                              if (ep_old.le.strain14) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(ep_old-strain13)

                          else
                              if (ep_old.le.strain15) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(ep_old-strain14)
                    else
                              if (ep_old.le.strain16) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(ep_old-strain15)
                    else
                              if (ep_old.le.strain17) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(ep_old-strain16)              
                    else
                              if (ep_old.le.strain18) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(ep_old-strain17)  
                    else
                              if (ep_old.le.strain19) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(strain18-strain17)
     9                       + stif193*(ep_old-strain18)
                    
                    else
                              if (ep_old.le.strain20) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(strain18-strain17)
     9                       + stif193*(strain19-strain18)
     1                       + stif203*(ep_old-strain19)
                  
                    else
                              if (ep_old.le.strain21) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(strain18-strain17)
     9                       + stif193*(strain19-strain18)
     1                       + stif203*(strain20-strain19)
     2                       + stif213*(ep_old-strain20)               
                    else
                              if (ep_old.le.strain22) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(strain18-strain17)
     9                       + stif193*(strain19-strain18)
     1                       + stif203*(strain20-strain19)
     2                       + stif213*(strain21-strain20) 
     3                       + stif223*(ep_old-strain21) 
                  
                    else
                              if (ep_old.le.strain23) then
                    yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(strain18-strain17)
     9                       + stif193*(strain19-strain18)
     1                       + stif203*(strain20-strain19)
     2                       + stif213*(strain21-strain20) 
     3                       + stif223*(strain22-strain21)
     4                       + stif233*(ep_old-strain22)  
                              else
                      yield3 = stressy3 + stif13*strain1 
     1                       + stif23*(strain2-strain1)
     2                       + stif33*(strain3-strain2)
     3                       + stif43*(strain4-strain3)
     4                       + stif53*(strain5-strain4)
     5                       + stif63*(strain6-strain5)
     6                       + stif73*(strain7-strain6)
     7                       + stif83*(strain8-strain7)
     8                       + stif93*(strain9-strain8)
     9                       + stif103*(strain10-strain9)
     1                       + stif113*(strain11-strain10)
     2                       + stif123*(strain12-strain11)
     3                       + stif133*(strain13-strain12)
     4                       + stif143*(strain14-strain13)
     5                       + stif153*(strain15-strain14)
     6                       + stif163*(strain16-strain15)
     7                       + stif173*(strain17-strain16)
     8                       + stif183*(strain18-strain17)
     9                       + stif193*(strain19-strain18)
     1                       + stif203*(strain20-strain19)
     2                       + stif213*(strain21-strain20) 
     3                       + stif223*(strain22-strain21)
     4                       + stif233*(strain23-strain22)
     5                       + stif233*(ep_old-strain23)      
                                        endif
                  endif
                    endif
                    endif
                    endif
                     endif
                      endif
                    endif
                    endif
                    endif
                    endif
                    endif
                 endif
                 
                                          endif
                                      endif
                                  endif
                              endif
                          endif
                      endif
                  endif              
              endif
          endif
          endif
          endif  
      r11 = yield_old/yield_old
      r22 = r11
      r33 = yield2/yield_old
      r12 = one/(4.d0/3.d0-one/3.d0*(yield_old/yield2)**two)**half
      r23 = one/(4.d0/3.d0*(yield_old/yield3)**two
     1      -one/3.d0)**half
      r13 = r23

C     Compute parameters F, G, H, L, M, N
      F = half*(one/r22**two+one/r33**two-one/r11**two)
      G = half*(one/r33**two+one/r11**two-one/r22**two)
      H = half*(one/r11**two+one/r22**two-one/r33**two)
      O = threehalfs/r23**two
      P = threehalfs/r13**two
      Q = threehalfs/r12**two
C Compute equivalent stress   
      yieldeff = (F*(stress22-stress33)**two+
     1            G*(stress33-stress11)**two+
     2            H*(stress11-stress22)**two+
     3            two*O*stress23**two+
     4            two*P*stress31**two+
     5            two*Q*stress12**two)**half
C     Yield function
      fyield = yieldeff - yield_old
              
      if (fyield .le. tol) then
C     Elastic range 
      stateNew(i,1)  = ep_old
      else
C     Plastic range
C     Flow direction
c      real*8 mn11, mn12, mn1y, mm1    
      z1 = (-G*(stress33 - stress11)+H*(stress11 - stress22))/yieldeff
      z2 =  (F*(stress22 - stress33)-H*(stress11 - stress22))/yieldeff
      z3 = (-F*(stress22 - stress33)+G*(stress33 - stress11))/yieldeff
      z4 =  two*Q*stress12/ yieldeff
      z5 =  two*O*stress23/ yieldeff
      z6 =  two*P*stress31/ yieldeff

      
      HF = z1**two*c11 + z2**two*c22 + z3**two*c33 
     1   + two*z1*z2*c12 + two*z2*z3*c23 + two*z1*z3*c13
     2   + z4**two*c44 + z5**two*c55 + z6**two*c66       
      
      
      delta_gamma = fyield/HF

C     Compute plastic modulus Hp
  
C First estimate delta_gamma
C Update stress
      stress11 = stress11 - delta_gamma*(c11*z1+c12*z2+c13*z3)
      stress22 = stress22 - delta_gamma*(c21*z1+c22*z2+c23*z3)
      stress33 = stress33 - delta_gamma*(c31*z1+c32*z2+c33*z3)
      stress12 = stress12 - delta_gamma*z4*c44
      stress23 = stress23 - delta_gamma*z5*c55
      stress31 = stress31 - delta_gamma*z6*c66

      

 
C Update plastic strain
      dstrain11 = delta_gamma * z1
      dstrain22 = delta_gamma * z2
      dstrain33 = delta_gamma * z3
      dstrain12 = delta_gamma * z4
      dstrain23 = delta_gamma * z5
      dstrain31 = delta_gamma * z6
      tr_dstrain = (dstrain11 + dstrain22 + dstrain33)/3.d0
      dev11 = dstrain11 - tr_dstrain
      dev22 = dstrain22 - tr_dstrain
      dev33 = dstrain33 - tr_dstrain
      dstrain   = sqrt23*(dev11**two+dev22**two+dev33**two
     1            + two*dstrain12**two+two*dstrain23**two
     2            + two*dstrain31**two)**half

      stateNew(i,1) = ep_old + dstrain
      stateNew(i,2) = delta_gamma

      
      endif
      stressNew(i,1) = stress11
      stressNew(i,2) = stress22
      stressNew(i,3) = stress33
      stressNew(i,4) = stress12
      stressNew(i,5) = stress23
      stressNew(i,6) = stress31     
C
C     Calculate stress state      
      sh123=(stress11+stress22+stress33)/3.0d0
      sig1=stress11-sh123
      sig2=stress22-sh123
      sig3=stress33-sh123
      smises1=(stress11-stress22)**2.d0+(stress22-stress33)**2.d0
     1        +(stress33-stress11)**2.d0
      smises2=(stress12**2.d0+stress23**2.d0+stress31**2.d0)
      smises3=0.5d0*(smises1+6.0d0*smises2)
      smises4=sqrt(smises3)
      smises =smises4+1.d0
      tri = sh123/(smises + 1.d0)
      J31=sig1*sig2*sig3+two*stress12*stress23*stress31
      J32=sig1*(stress23**two)+sig2*(stress31**two)+sig3*(stress12**two)
      J3=J31-J32
      xi1=J3/(smises**3.d0+1.d0)
      xi=27.d0*xi1/two
      stateNew(i,3) = tri
      stateNew(i,4) = xi

C Update the specific internal energy -
        enerInternNew(i) = enerInternOld(i)
C Update the dissipated inelastic specific energy -
        enerInelasNew(i) = enerInelasOld(i)
  100 continue
C
      return
      end