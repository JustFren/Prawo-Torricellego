import math

def euler_method(x0,y0,epsilon,x,func):
    #x0 is the starting point,
    #y0 is the initial value
    #x is the target 
    #epsilon is the precision
    #func is the right side of a first order ODE
    #ex.:
    #   dx/dt = 2x
    y=y0
    while x0<x:
        temp=y
        y=y+epsilon*func(x0,y)
        x0+=epsilon

    return y

def func(x,y):
    return y

print(euler_method(0,1,0.00001,1,func))