import gepard as g
from gepard.fits import th_KM10b
import math as m
import random
import sys
import warnings

def main(args):
    xsec_mean=0
    phi_mean=0
    xB_mean=0
    t_mean=0
    Q2_mean=0
    N=1000
    n=0

    while(n<N):
        xB_val=random.uniform(args[4], args[5])
        t_val=random.uniform(max(args[0],-1.0), args[1])
        Q2_val=random.uniform(args[2], min(args[3],2.0*0.938*xB_val*11))
        phi_val=random.uniform(args[6], args[7])*m.pi/180.

        pt = g.DataPoint(xB=xB_val, t=t_val, Q2=Q2_val, phi=phi_val, process="ep2epgamma", exptype='fixed target', in1energy=10.6, in1charge=-1,in1polarization=0,frame='Trento')
        pt.to_conventions()
        xsec_val=th_KM10b.XS(pt)

        if(not m.isnan(xsec_val)):
            n+=1
            xsec_mean+=1.0/xsec_val
            #print(phi_rad, xsecval)

    xB_mean=args[8] #xB_mean/N
    t_mean=args[9] #t_mean/N
    Q2_mean=args[10] #Q2_mean/N
    phi_mean=args[11]*m.pi/180. #phi_mean/N

    xsec_mean=n/xsec_mean
    pt = g.DataPoint(xB=xB_mean, t=t_mean, Q2=Q2_mean, phi=phi_mean, process="ep2epgamma", exptype='fixed target', in1energy=10.6, in1charge=-1,in1polarization=0,frame='Trento')
    pt.to_conventions()
    xsec_point=th_KM10b.XS(pt)
    
    #print(xsec_mean, xsec_point,phi_mean,xB_mean,t_mean,Q2_mean)
    return xsec_mean/xsec_point

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    args = [float(x) for x in sys.argv[1:]]     
    x =main(args)
    print(x)
