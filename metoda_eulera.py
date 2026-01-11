import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import CubicSpline
from approximations_of_a import approximated_a_10

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


times=[0,64_800,151_200,241_200,327_600,414_000,583_200,669_600,756_000,842_400,928_800,1_008_000]
volumes=[147,134,115,102,85,63,44,24,17,12,8,4]
def compare_methods(x0,x,y0,epsilon,func,draw=1):
    euler=euler_method_array(x0,x,epsilon,y0,func)
    t_space = np.linspace(x0, x, len(euler))
    lsode=scipy.integrate.odeint(func,y0,t_space,tfirst=True)
    
    if draw==0:
        return [euler,lsode]
    err=calc_error(euler,lsode)
    obj=[147,134,115,102,85,63,44,24,17,12,8,4]
    for i in range(len(obj)):
        obj[i]=volume_to_height(obj[i]*(10**6))


    x=np.linspace(x0,x,len(euler))
    fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,sharex=True)
    #ax.set_yscale("log")
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("Water level (m)")
    ax1.plot(t_space/3600,euler,label="Euler", linewidth = 3.5)
    ax1.plot(t_space/3600,lsode,label="ODEINT/LSODE", color = "orange", linewidth=1)
    
    ax1.plot(np.array(times)/3600,obj,"^",label="Real data", color = "red")
    ax1.legend()
    ax2.set_xlabel("Time (h)")
    ax2.set_ylabel("Error between Euler and LSODE (Euler-LSODE)")
    ax2.plot(t_space/3600,err)


    #print(f"average error = {average_err}\nmean error = {mean_err}")
    plt.show()
    return [euler,lsode]

def compare_methods_mult(times,volumes,epsilon,func):
    volumes=[volume_to_height(volumes[i]*(10**6)) for i in range(len(volumes))]
    euler=[]
    lsode=[]
    for i in range(len(times)-1):
        euler.append(euler_method_array(times[i],times[i+1],epsilon,volumes[i],func))
        lsode.append(scipy.integrate.odeint(func,volumes[i],np.linspace(times[i],times[i+1],len(euler)),tfirst=True))
    
    #err=calc_error(euler,lsode)
    obj=[147,134,115,102,85,63,44,24,17,12,8,4]
    for i in range(len(obj)):
        obj[i]=volume_to_height(obj[i]*(10**6))


    fig, ax1 = plt.subplots(nrows=1,ncols=1,sharex=True)
    #ax.set_yscale("log")
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("Water level (m)")
    #ax2.set_xlabel("Time (s)")
    #ax2.set_ylabel("Error between Euler and LSODE (Euler-LSODE)")
    for i in range(len(times)-1):
        x=np.linspace(times[i],times[i+1],len(euler[i]))
        ax1.plot(x/3600,euler[i],label="Euler",color="blue")
        #ax1.plot(x/3600,lsode[0][i],label="ODEINT/LSODE")
    
    ax1.plot(np.array(times) / 3600,obj,"^",label="Real data",color="red")
    #ax1.legend()

    #ax2.plot(x,err)
    


    #print(f"average error = {average_err}\nmean error = {mean_err}")
    plt.show()
    #return [euler,lsode]

#dla wody mi pomiędzy [0.60,0.65]
mi=0.62
a=612 # 648
A0=26_300_000
g=9.81
h0=24.41
Htotal=35.41
def A(h):
    return A0 * (h/Htotal)**2
Vwirt=1/3*A(h0)*h0

def height_to_volume(h):
    h=h+h0
    return ((A(h)*(h)/3)-Vwirt)

def volume_to_height(V):
    return np.cbrt(3*(V+Vwirt)*(Htotal**2)/A0)-h0

inflow=[1040,644.1,408.2,302.5,230.8,172.5,130.5,116.5,104.5,92.2,81.3,82.8]
wsp=mi*a
G=2*g
Q_in=CubicSpline(times,inflow)


apprx_a_func=CubicSpline(times,approximated_a_10)
def approx_raciborz(t,h):
    def Q_out(t,h):
        return mi*a*apprx_a_func(t)*np.sqrt(G*h)
    if h<=0:
        return 0
    return (Q_in(t)-Q_out(t,h))/A(h+h0)

compare_methods(0,times[-1],volume_to_height(volumes[0]*(10**6)),1,approx_raciborz)  
compare_methods_mult(times,volumes,1,approx_raciborz)

#wynik=compare_methods(0,1_000_000,volume_to_height(147_000_000),1,approx_raciborz)

testlist=[147,134,115,102,85,63,44,24,17,12,8,4]

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

def approx_raciborz_redundant(t,h):
    def Q_out(h):
        return min(1200,wsp*np.sqrt(G*h))  
    if h<=0:
        return 0
    return (Q_in(t)-Q_out(h))/A(h+h0)


def model_prognostyczny(times, volumes_imgw, q_in_list, epsilon=1):

    V_max = 185000000
    V_80 = 0.8 * V_max
    Q_safe = 1210.0
    
    q_in_func = CubicSpline(times, q_in_list)
    

    t_start = times[0]
    t_end = times[-1]
    total_seconds = t_end - t_start
    
    t_sim = np.linspace(t_start, t_end, total_seconds + 1)
    v_sim = np.zeros(len(t_sim))
    
    v_sim[0] = volumes_imgw[0] * 1e6 
    
    for i in range(total_seconds):
        t_teraz = t_sim[i]
        v_teraz = v_sim[i]
        q_in = q_in_func(t_teraz)
        
        h_teraz = volume_to_height(v_teraz)
        if h_teraz > 0:
            q_max = wsp * np.sqrt(2 * g * h_teraz)
        else:
            q_max = 0

        if v_teraz < V_80:
            if v_teraz > 0:
                q_out = Q_safe
            else:
                q_out = q_in
        else:
            q_out = max(Q_safe, q_in)
        
        q_out = min(q_out, q_max)
            
        dv = (q_in - q_out) * epsilon
        v_sim[i+1] = max(0, v_teraz + dv)
    

    plt.plot(t_sim / 3600, v_sim / 1e6, color='blue', linewidth=2, label="Prognoza modelu")
    
    plt.plot(np.array(times) / 3600, volumes_imgw, "^", color='red', label="Dane rzeczywiste IMGW", markersize=10)


    plt.axhline(y=148, color='orange', linestyle=':', label="Próg 80% (148 mln m³)")
    plt.fill_between(t_sim / 3600, 0, v_sim / 1e6, color='blue', alpha=0.1)

    plt.title("Symulacja objętości zbiornika", fontsize=14)
    plt.xlabel("Czas od startu [h]", fontsize=12)
    plt.ylabel("Objętość [$mln/m^3$]", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()

model_prognostyczny(times, volumes, inflow, 1)