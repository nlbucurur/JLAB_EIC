import gepard as g
from gepard.fits import th_KM15
import math as m
import random
import sys
import warnings

def main(args):
    N=args[4]
    phi_val=m.pi - args[6]*m.pi/180.
    
    pt = g.DataPoint(xB=args[0], Q2=args[1], t=args[2], phi=phi_val, process="ep2epgamma", exptype='fixed target', in1energy=10.6, in1charge=-1,in1polarization=0,in2polarization='U',frame='Trento')
    threed=th_KM15.XS(pt)

    pt = g.DataPoint(xB=args[3], Q2=args[4], t=args[5], phi=phi_val, process="ep2epgamma", exptype='fixed target', in1energy=10.6, in1charge=-1,in1polarization=0,in2polarization='U',frame='Trento')
    fourd=th_KM15.XS(pt) 
    
    return fourd/threed

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    args = [float(x) for x in sys.argv[1:]]     
    x =main(args)
    print(x)
