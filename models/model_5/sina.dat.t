   0   0   0   0   0   0   0  = ifvolt,ifadat,ifaeg,ifb0,iefrr,ifmove,ifwoff
  200   0   0   0   0   0   0  = nstep,nwave,ifconv,iftime,ifcont,ifbdat,ifje
 1.00000E-06 0.00000E+00   0  = dt(s), time0(s), istep2
 0.00000E+00 0.00000E+00      = freq(Hz), freq-end(Hz)
   0  15   0                  = ifmat, mite, ia18
   0
 1.00000E-03 1.00000E-02      = berr(T), berr2(T)
 1.00000E+00                  = 加速係数
 1.00000E-06 0.00000E+00      = ddmin(収束判定値), ddmin(閾値)(%)
   1                          = ntnod(測定節点数)
  1              = 節点番号(6i8)
   0                          = nrotmat                 ! 回転材質数
                              = rotmat(1),rotmat(2).... ! 材質番号
 0.00000E+00 2.09440E+00 4.18880E+00 1.00000E+04     = VMIN,VMAX,VMIDOL,VA
 1.50000E+01 0.00000E+00 3.60000E+01                 = PERIOD, zure, test
   0   1                                             = ifabs, ia18@femeeme
   0


*ia18