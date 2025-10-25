# L.R.M./M.B./S.G./A.H. 07/2021
# ORCA V. 4.2.1
# DLPNO-CCSD(T)/mixedCBS//B97-3c or DLPNO-CCSD(T1)/CBS(T|vt/Q|vt)//B97-3c
# T = def2-TZVPP
# Q = def2-QZVPP
# CBS extrapolation according to Neese/Valeev
# t = TightPNO settings
# vt = VeryTightPNO settings
# mixedCBS = CBS(T|vt/Q|t+DS)+DT1
# T(VeryTightPNO)/Q(TightPNO) extrapolation
# DS = Ecorr(DLPNO-CCSD(T)/def2-SVP|vt)-Ecorr(DLPNO-CCSD(T)/def2-SVP|t)
# DT1 = Ecorr(DLPNO-CCSD(T1)/def2-SVP|t)-Ecorr(DLPNO-CCSD(T)/def2-SVP|t)
if [ "$TMER" == "" ]; then
  tmer=tmer2++
else
  tmer=$TMER
fi

f=$1
w=$2
$tmer m1/$f m2/$f m3/$f						x -1 -1  1		$w -40.05   1
$tmer m3/$f m4/$f						x -1  1			$w -37.71   2
$tmer m1/$f m5/$f m4/$f m6/$f					x -1 -1  1  1		$w -12.59   3
$tmer m7/$f m8/$f m9/$f						x -1 -1  1		$w -17.04   4
$tmer m9/$f m8/$f m10/$f					x -1 -1  1		$w -10.66   5
$tmer m10/$f m11/$f						x -1  1			$w -12.66   6
$tmer m11/$f m8/$f m12/$f					x -1 -1  1		$w -14.73   7
$tmer m12/$f m13/$f						x -1  1			$w -10.51   8
$tmer m14/$f m13/$f						x -1  1			$w -12.76   9
$tmer m14/$f m15/$f						x -1  1			$w -21.20   10
$tmer m13/$f m15/$f						x -1  1			$w -8.44    11
$tmer m9/$f m16/$f m15/$f m8/$f 				x -1 -1  1  1		$w -2.79    12
$tmer m19/$f m20/$f m21/$f m22/$f m23/$f			x -1 -1 -1  1  1	$w -178.96  13
$tmer m24/$f m25/$f m26/$f					x -1 -1  1		$w -42.33   14
$tmer m26/$f m27/$f m28/$f m29/$f				x -1 -1  1  1		$w -7.17    15
$tmer m30/$f m31/$f m32/$f					x -1 -1  1		$w -193.43  16
$tmer m33/$f m31/$f m34/$f					x -1 -1  1		$w -200.10  17
$tmer m35/$f m31/$f m36/$f					x -1 -1  1		$w -203.09  18
$tmer m37/$f m38/$f m39/$f					x -1 -0.5 1		$w -4.94    19
$tmer m40/$f m41/$f m42/$f					x -1 -1  1		$w -27.02   20
$tmer m43/$f m44/$f m45/$f m46/$f				x -1 -1  1  1		$w -25.09   21
$tmer m43/$f m46/$f m47/$f					x -1 -1  1		$w -46.70   22
$tmer m49/$f m50/$f m51/$f					x -1 -1  1		$w -115.35  23
$tmer m51/$f m52/$f m53/$f					x -1 -1  1		$w -46.37   24
$tmer m54/$f m55/$f m56/$f m17/$f m41/$f			x -1 -2  1  2  2	$w -19.39   25
$tmer m56/$f m57/$f m58/$f m18/$f m59/$f 			x -1 -1  1  1  1	$w -47.73   26
$tmer m60/$f m48/$f m61/$f m62/$f				x -1 -1  1  1		$w -25.86   27
$tmer m63/$f m48/$f m64/$f m62/$f				x -1 -1  1  1		$w -28.77   28
$tmer m65/$f m66/$f						x -1  1			$w -11.13   29
$tmer m67/$f m41/$f m65/$f					x -1 -1  1		$w -141.32  30
$tmer m68/$f m31/$f m69/$f					x -1 -1  1		$w -63.73   31
$tmer m70/$f m20/$f m71/$f m41/$f				x -1 -1  1  1		$w -49.75   32
$tmer m72/$f m73/$f m74/$f					x -1  1  1		$w -10.64   33
$tmer m75/$f m76/$f m77/$f m17/$f				x -1 -1  1  1		$w -0.65    34
$tmer m78/$f m77/$f						x -1  1			$w -5.46    35
$tmer m79/$f m80/$f m81/$f m21/$f				x -1 -4  1  4		$w -49.04   36
$tmer m82/$f m83/$f m84/$f					x -1 -1  1		$w -40.62   37
$tmer m85/$f m86/$f m87/$f					x -1 -1  1		$w -27.70   38
$tmer m88/$f m89/$f m90/$f m91/$f				x -1 -1  1  1		$w -44.99   39
$tmer m92/$f m93/$f m94/$f m91/$f				x -1 -1  1  1		$w -24.51   40
$tmer m94/$f m95/$f m96/$f m93/$f				x -1 -1  1  1		$w -5.96    41
$tmer m96/$f m97/$f m98/$f m99/$f				x -1 -1  1  0.5		$w -26.96   42
$tmer m100/$f m101/$f m102/$f m17/$f m103/$f m104/$f m105/$f	x -1 -1 -1 -2  1  1  1	$w -80.58   43
$tmer m106/$f m102/$f m107/$f m101/$f m108/$f m104/$f m105/$f	x -1 -1 -1 -1  1  1  1	$w -66.15   44
$tmer m109/$f m83/$f m110/$f					x -1 -1  1		$w -42.83   45
$tmer m111/$f m112/$f m113/$f					x -1 -1  1		$w -33.27   46
$tmer m114/$f m115/$f m116/$f m41/$f m107/$f			x -1 -1  1  2  1	$w -58.47   47
$tmer m117/$f m91/$f m118/$f					x -1 -1  1		$w -4.24    48
$tmer m119/$f m120/$f m121/$f m41/$f				x -1 -1  1  1		$w -35.26   49
$tmer m121/$f m122/$f m91/$f					x -1  1  1		$w -2.74    50
$tmer m122/$f m123/$f m124/$f					x -1 -1  1		$w -29.77   51
$tmer m125/$f m126/$f						x -1  1			$w -42.72   52
$tmer m127/$f m120/$f m128/$f m41/$f				x -1 -1  1  1		$w -30.78   53
$tmer m129/$f m130/$f m131/$f					x -1 -1  1		$w -12.33   54
$tmer m132/$f m133/$f						x -1  1			$w -4.61    55
$tmer m134/$f m135/$f m136/$f m137/$f				x -1 -1  1  1		$w -13.65   56
$tmer m138/$f m139/$f						x -1  1			$w -31.89   57
$tmer m140/$f m141/$f						x -1  1			$w -39.30   58
$tmer m142/$f m143/$f m144/$f					x -1 -1  1		$w -29.89   59
$tmer m145/$f m146/$f m147/$f					x -1 -1  1		$w -66.00   60
$tmer m148/$f m149/$f m150/$f					x -1 -1  1		$w -69.48   61
