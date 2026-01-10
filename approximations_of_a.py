from scipy.interpolate import CubicSpline
import numpy as np
#dla wody mi pomiędzy [0.60,0.65]
mi=0.62
a=612 # 648
A0=26_300_000
g=9.81
h0=24.41
Htotal=35.41
inflow=[1040,644.1,408.2,302.5,230.8,172.5,130.5,116.5,104.5,92.2,81.3,82.8]
times=[0,64_800,151_200,241_200,327_600,414_000,583_200,669_600,756_000,842_400,928_800,1_008_000]
volumes=[147,134,115,102,85,63,44,24,17,12,8,4]
wsp=mi*a/2
G=2*g
Q_in=CubicSpline(times,inflow)

def A(h):
    return A0 * (h/Htotal)**2
def _volume_to_height(V):
    Vwirt=1/3*A(h0)*h0
    return np.cbrt(3*(V+Vwirt)*(Htotal**2)/A0)-h0

def approx_a(i,prec):
    approx_a_mult=0.5
    #approx_a_mult is a value between 0 and 1 and it indicates how much of the outflow capacity is being used
    def _approx_raciborz(t,h,approx_a_mult):
        wsp=mi*a*approx_a_mult 
        if h<=0:
            return 0
        return (Q_in(t)-wsp*np.sqrt(G*h))/A(h+h0)
    def _euler_method(x0,x,y0,epsilon,func,approx_a_mult):
        y=y0
        while x0<x:
            y=y+epsilon*func(x0,y,approx_a_mult)
            x0+=epsilon
        return y
    end_val=0
    temp_volumes=[_volume_to_height(volumes[i]*(10**6)) for i in range(len(volumes))]
    for j in range(prec):
        end_val=_euler_method(times[i],times[i+1],temp_volumes[i],1,_approx_raciborz,approx_a_mult)
        #print(end_val)
        if end_val>temp_volumes[i+1]:
            approx_a_mult+=(1/2)**(j+2)
        elif end_val<temp_volumes[i+1]:
            approx_a_mult-=(1/2)**(j+2)
    return approx_a_mult

approximated_a_10=[0.21142578125,0.21142578125,0.15771484375,\
                0.11376953125,0.11279296875,0.12255859375,0.08056640625,0.13623046875,\
                0.09228515625,0.08837890625,0.09033203125,0.11767578125]

if __name__ == "__main__":
    for i in range(len(times)-1):
        print(approx_a(i,10))