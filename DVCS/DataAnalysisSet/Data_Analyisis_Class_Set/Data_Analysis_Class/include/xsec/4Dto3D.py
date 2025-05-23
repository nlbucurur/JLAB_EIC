import gepard as g
from gepard.fits import th_KM15
import math as m
import random
import sys
import warnings

def main(args):    
    pt = g.DataPoint(xB=args[0], Q2=args[1], t=args[2], phi=args[3]*m.pi/180., process="ep2epgamma", exptype='fixed target', in1energy=10.6, in1charge=-1,in1polarization=0,in2polarization='U',frame='Trento')
    pt.to_conventions()
    threed=th_KM15.XS(pt)

    pt = g.DataPoint(xB=args[4], Q2=args[5], t=args[6], phi=args[7]*m.pi/180., process="ep2epgamma", exptype='fixed target', in1energy=10.6, in1charge=-1,in1polarization=0,in2polarization='U',frame='Trento')
    pt.to_conventions()
    fourd=th_KM15.XS(pt) 
    
    return fourd/threed

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    args = [float(x) for x in sys.argv[1:]]     
    x =main(args)
    print(x)
