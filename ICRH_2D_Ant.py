import math as m
import SI
cm = 10.0**-2
PI = m.pi

def magneticField(x,B0,R0):
    return (B0*R0)/x

def def_params():
    L_plasm = 35.0*cm # 35*cm %19*cm
    N0 = 1.0 * 10.0**18 # % m^-3 L-mode: N0 = 1.0 * 10.^18; dx = 4.8*cm;
    dx = 1.8*cm  #1*cm
    w = 11.8*cm  #
    wz = w/2     #
    s = 1.6*cm   #
    d = 3.8*cm   #
    t = 0.94*cm  #
    return N0,dx,w,wz,s,d,t,L_plasm


def plasmParam_GlobusM1():
    R0 = 0.36
    a  = 0.24
    B0 = 0.4
    R_lfs_bnd = R0 + a
    B  = magneticField(R_lfs_bnd,B0,R0)
    w_ci = (SI.e0*B)/(2*SI.m_PROTO)
    return w_ci

def plasmParam_GlobusM2():
    R0 = 0.36
    a  = 0.24
    B0 = 1.0
    R_lfs_bnd = R0 + a
    B  = magneticField(R_lfs_bnd,B0,R0)
    w_ci = (SI.e0*B)/(2*SI.m_PROTO)
    return w_ci

def alpha(f,N0,dx, w_ci):
    α = (1.95*10**-6)*( 0.5*(N0/dx)*(2*PI*f/w_ci)**2 )**(1/3)
    return α
def delta(alpha,dx,d):
    Δ = 1 + alpha*(dx+d) + (alpha*(dx+d))**2
    return Δ

def N_(alpha, dx, d):
    res = d*alpha*( 0.5 + alpha*(dx+d) )
    return res
def f_q(m,wz,t):
    res = 0.5*m*PI/(wz + t);
    return res




# # END GAME CALCULATION: series_C series_L series_R
def AntennaICRH__ser_C(s,d,wz,t):
    Sum=0.0
    for j in range(0,3): # 0:2
        pm = j*2+1
        q = f_q(pm,wz,t)
        Sum=Sum - (1/SI.EPS0/pm/PI) *( (1-m.exp(2*q*s))*(1-m.exp(-2*q*d))/(m.exp(2*q*s)-m.exp(-2*q*d)) )* (m.sin(q*wz)/q/wz)**2
    C = 1.0/Sum
    return C
def AntennaICRH__ser_L(wz,d,N,del11):
    L = 0.5*SI.MU0*d/wz*( 1 - N/del11 )
    return L
def AntennaICRH__ser_R(f, wz, d, alpha, del11):
    R = (0.5*SI.MU0*2*PI*f/wz)*(0.87*d**2.0*alpha/del11)
    return R

