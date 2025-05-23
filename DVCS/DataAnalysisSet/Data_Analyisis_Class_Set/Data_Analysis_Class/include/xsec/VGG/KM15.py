from gepard.cff import (CFF, MellinBarnesCFF, DispersionCFF, PionPole,  # noqa: F401
                  DispersionFixedPoleCFF, DispersionFreePoleCFF,
                  HybridFixedPoleCFF, HybridFreePoleCFF, GoloskokovKrollCFF)
from gepard.data import (DataPoint, DataSet, dset, loaddata,  # noqa: F401
                   select, list_data, describe_data)
from gepard.dis import DIS  # noqa: F401
from gepard.dvcs import DVCS, BM10, BMK, BM10ex, BM10tw2, hotfixedBMK  # noqa: F401
from gepard.dvmp import DVMP, MellinBarnesTFF  # noqa: F401
from gepard.eff import ZeroEFF, DipoleEFF, KellyEFF  # noqa: F401
from gepard.fitter import MinuitFitter  # noqa: F401
from gepard.gpd import GPD, ConformalSpaceGPD, TestGPD, PWNormGPD  # noqa: F401
from gepard.kinematics import tmin, tmax, weight_BH, prepare  # noqa: F401
from gepard.qcd import beta, as2pf  # noqa: F401
from gepard.theory import Theory  # noqa: F401
from gepard.fits import th_KM15

import sys

def main(args):
    pt = DataPoint(xB=0.30441159813051266, t=-0.339796, Q2=5.04437,Qp2=0.)
    pt.to_conventions()    
    test=th_KM15.ImH(pt),
    
    pt = DataPoint(xB=args[0], t=args[1], Q2=args[2], Qp2=args[3])
    pt.to_conventions()
    return [th_KM15.ReH(pt), th_KM15.ImH(pt), th_KM15.ReE(pt), th_KM15.ImE(pt), th_KM15.ReHt(pt), th_KM15.ImHt(pt), th_KM15.ReEt(pt), th_KM15.ImEt(pt)]  # This will work
    #pt.xi = 0.0135
    #pt.xip= -0.012
    #print(pt.xi, pt.xip)
    #return [th_KM15.ReH(pt), th_KM15.ImH(pt), th_KM15.ReE(pt), th_KM15.ImE(pt), th_KM15.ReHt(pt), th_KM15.ImHt(pt), th_KM15.ReEt(pt), th_KM15.ImEt(pt)]  # This will work
    #return [th_KM15.cff(pt)[0], th_KM15.cff(pt)[1], th_KM15.cff(pt)[2], th_KM15.cff(pt)[3], 0, 0, 0, 0]  # This will work
    

if __name__ == "__main__":
    args = [float(x) for x in sys.argv[1:]]     
    x =main(args)
    print(x)

