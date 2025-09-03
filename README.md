# Anisotropic-VUMAT

Required material:
*DENSITY
*DEPVAR = 2
*USER MATERIAL
  props(1) = Elastic modulus of coupon at 0 degree
  props(2) = Elastic modulus of coupon at 90 degree
  props(3) = In-plane Poisson's ratio 
  props(4) = Poisson's ratio measured coupons at 90 degree
  props(5) = In-plane shear modulus
  props(6) = Normal-plane shear modulus
  
  props(7) = True yield stress of coupon at 0 degree
  props(8-30) = Strain hardening rate between predifined true plastic strain1 to strain23 of coupon at 0 degree
  
  props(31) = True yield stress of coupon at 90 degree
  props(32-54) = Strain hardening rate between predifined true plastic strain1 to strain23 of coupon at 90 degree
 
  props(55) = True yield stress of coupon at 45 degree
  props(56-78) = Strain hardening rate between predifined true plastic strain1 to strain23 of coupon at 45 degree
