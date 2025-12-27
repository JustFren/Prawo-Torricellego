import math

def euler_method(x0,y0,epsilon,x,func):
    #x0 is the starting point,
    #y0 is the initial value
    #x is the target 
    #epsilon is the precision
    #func is the right side of a first order ODE
    #ex.:
    #   dy/dx = 2y => y=e^2x
    y=y0
    while x0<x:
        temp=y
        y=y+epsilon*func(x0,y)
        x0+=epsilon

    return y

def func(x,y):
    return y
#przybliżenie liczby eulera
#print(euler_method(0,1,0.00001,1,func))

#======================================================
#prawo Torricellego dla cylindra

#promień cylindra
r=1
#pole powierzchni przekroju cylindra
A=math.pi*(r**2)
#pole powierzchni otworu wylotowego na spodzie cylindra
a=0.01*r
#przyspieszenie grawitacyjne
g=9.81

def T_cylinder(t,h):
    if h<0:
        return 0
    return -(a/A)*math.sqrt(2*g*h)

#przybliżenie poziomu zbiornika w czasie t w zależności od poziomu startowego h0
h0=100
t0=0
t=3600
print(euler_method(t0,h0,0.01,t,T_cylinder))
#=======================================================
#prawo Torricellego dla stożka/lejka

#promień stożka
R0=100
R1=10
#wysokość stożka
H=30
#przyspieszenie grawitacyjne
g=9.81

# -v- wzięte z Eur. J. Phys. 42 (2021) 065808 (11pp) strona 4
def T_funnel(t,h):
    if h<0:
        return 0
    return -math.sqrt(2*g*h)/((1+((R1-R0)*h)/(R0*H))**2)
t0=0
h0=H
t=100

print(euler_method(t0,h0,0.0001,t,T_funnel))
#TODO jakie przybliżenie zbiornika wybrać