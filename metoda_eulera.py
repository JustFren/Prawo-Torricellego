import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
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
    average_err=sum(err)/len(err)
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
    ax2.set_ylabel("Error between Euler and LSODE")
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
Vwirt=A0*(h0**3)/(3*(Htotal**2))
def A(h):
    return A0 * (h/Htotal)**2

wsp=(-1)*mi*a
G=2*g
def approx_raciborz(t,y):
    if y<=0:
        return 0
    return (wsp * np.sqrt(G*y))/A(y+h0)


def volume_to_height(V):
    return math.cbrt(3*(V+Vwirt)*(Htotal**2)/A0)-24.41

testlist=[147_000_000,134_000_000,128_000_000,115_000_000,57_000_000,21_000_000,17_000_000,12_000_000,10_000_000,0]

for i in range(len(testlist)):
    print(volume_to_height(testlist[i]))




wynik=compare_methods(0,64_000,volume_to_height(147_000_000),1,approx_raciborz)

for i in range(len(testlist)):
    testlist[i]=volume_to_height(testlist[i])
plot(testlist)
"""
* 17 września o 16:00 poziom wody wynosił ok. 147 mln m³,
* 18 września o 10.00 poziom wody wynosił ok. 134 mln m³,
* 19 września o 10.00 poziom wody wynosił już ok. 115mln m³ (ok. 65% pojemności zbiornika),
* 20 września o 11:00 poziom wody wynosił ok. 102 mln m³,
* 21 września o 11:00 poziom wody wynosił ok. 85 mln m³,
* 22 września o 11:00 poziom wody wynosił ok. 63 mln m³,
* 23 września o 10:00 poziom wody wynosił ok. 44 mln m³,
* 24 września o 10:00 poziom wody wynosił ok. 24 mln m³,
* 25 września o 10:00 poziom wody wynosił ok. 17 mln m³,
* 26 września o 10:00 poziom wody wynosił ok. 12 mln m³,
* 27 września o 10:00 poziom wody wynosił ok. 8 mln m³,
* 28 września o 10:00 poziom wody wynosił ok. 4 mln m³,
* 29 września 2024 ok godziny 8 zakończyła się praca zbiornika Racibórz Dolny, poziom wody wynosi 0m³.
"""
#TODO zmienić objętości na wartości h
