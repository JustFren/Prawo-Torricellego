import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import CubicSpline
def plot(list):
    # make data
    #x = np.linspace(x0, x, int((x-x0)/epsilon))
    x = [i for i in range(0,len(list))]
    y=list
    # plot
    fig, ax = plt.subplots()
    ax.plot(x,y)
    plt.show()

def euler_method_array(x0,x,epsilon,y0,func,draw=0):
    #x0 is the starting point,
    #y0 is the initial value
    #x is the target 
    #epsilon is the precision
    #func is the right side of a first order ODE
    #ex.:
    #   dy/dx = 2y => y=e^2x
    y=y0
    out=[]
    out.append(y0)
    x0+=2*epsilon
    while x0<x:
        temp=y
        y=y+epsilon*func(x0,y)
        #print(epsilon*func(x0,y))
        out.append(y)
        x0+=epsilon

    if draw==1:
        plot(out,x0,x,epsilon)
    return out


def calc_error(list1,list2):
    out=[]
    for i in range(len(list1)):
        out.append(list1[i]-list2[i])
    return out

def compare_methods(x0,x,y0,epsilon,func):
    euler=euler_method_array(x0,x,epsilon,y0,func)
    lsode=scipy.integrate.odeint(func,y0,np.linspace(x0,x,len(euler)),tfirst=True)
    err=calc_error(euler,lsode)
    #mean_err=err
    #mean_err.sort()
    #Jmean_err=mean_err[len(mean_err)//2]


    x=np.linspace(x0,x,len(euler))
    fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,sharex=True)
    #ax.set_yscale("log")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Water level (m)")
    ax1.plot(x,euler,label="Euler")
    ax1.plot(x,lsode,label="ODEINT/LSODE")
    ax1.legend()
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Error between Euler and LSODE (Euler-LSODE)")
    ax2.plot(x,err)


    #print(f"average error = {average_err}\nmean error = {mean_err}")
    plt.show()
    return [euler,lsode]
#promień cylindra
r=10
#pole powierzchni przekroju cylindra
A=math.pi*(r**2)
#pole powierzchni otworu wylotowego na spodzie cylindra
a=0.06*A
#przyspieszenie grawitacyjne
g=9.81

def test_func(t,y):
    return y

#dla wody mi pomiędzy [0.60,0.65]
mi=0.62
a=612 # 648
A0=26_300_000
g=9.81
h0=24.41
Htotal=35.41
#======
#model nie jest do końca dokładny więc trzeba wprowadzić zmienną "terenową"
#według tego modelu maksymalna pojemność to 208 mln m^3, ta zmienna to poprawia
zm_ter=0.889
#====
def A(h):
    return A0 * (h/Htotal)**2
Vwirt=1/3*A(h0)*h0

def height_to_volume(h):
    return ((A(h)*(h)/3)-Vwirt)*zm_ter

def volume_to_height(V):
    return np.cbrt(3*(V+Vwirt)*(Htotal**2)/A0)

print(volume_to_height(134_000_000))



#========v
times=[0,64_800,151_200,241_200,327_600,414_000,583_200,669_600,756_000,842_400,928_800,1_008_000]
inflow=[1040,644.1,408.2,302.5,230.8,172.5,130.5,116.5,104.5,92.2,81.3,82.8]
cs=CubicSpline(times,inflow)


def Q_in(t,h):
    def Q_calc(h,inflow):
        return volume_to_height(inflow+height_to_volume(h))-volume_to_height(height_to_volume(h))
    #xs=np.arange(0,1_000_000,1)
    #plot(cs(xs))
    return Q_calc(h,cs(t))
#========^    

wsp=(-1)*mi*a*(1/2)
G=2*g
def approx_raciborz(t,h):
    if h<=0:
        return 0
    return Q_in(t,h)+((wsp * np.sqrt(G*h))/A(h+h0))

testlist=[147,134,115,102,85,63,44,24,17,12,8,4]

#for i in range(len(testlist)):
#    print(volume_to_height(testlist[i]))




wynik=compare_methods(0,64_800,volume_to_height(147_000_000),0.5,approx_raciborz)




for i in range(len(testlist)):
    testlist[i]=volume_to_height(testlist[i]*(10**6))
plot(testlist)
"""
* 17 września o 16:00 poziom wody wynosił ok. 147 mln m³,t=0s
* 18 września o 10.00 poziom wody wynosił ok. 134 mln m³,t=64_800s
* 19 września o 10.00 poziom wody wynosił już ok. 115mln m³ (ok. 65% pojemności zbiornika),t=151_200s
* 20 września o 11:00 poziom wody wynosił ok. 102 mln m³,t=241_200
* 21 września o 11:00 poziom wody wynosił ok. 85 mln m³,t=327_600
* 22 września o 11:00 poziom wody wynosił ok. 63 mln m³,t=414_000
* 23 września o 10:00 poziom wody wynosił ok. 44 mln m³,t=496_800
* 24 września o 10:00 poziom wody wynosił ok. 24 mln m³,t=583_200
* 25 września o 10:00 poziom wody wynosił ok. 17 mln m³,t=669_600
* 26 września o 10:00 poziom wody wynosił ok. 12 mln m³,t=756_000
* 27 września o 10:00 poziom wody wynosił ok. 8 mln m³,t=842_400
* 28 września o 10:00 poziom wody wynosił ok. 4 mln m³,t=928_800
* 29 września 2024 ok godziny 8 zakończyła się praca zbiornika Racibórz Dolny, poziom wody wynosi 0m³.t=1_008_000
"""
#TODO zmienić objętości na wartości h
