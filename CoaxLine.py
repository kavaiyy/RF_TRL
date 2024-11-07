import cmath as cm
import math as m
import SI

PI = cm.pi
# def coth(x):
#     return 1/cm.tanh(x)

# L_coax in F/m
# input in SI
# C_coax in F/m
# input is SI  

def serCoax_L(r,R):
    L = (SI.MU0/2.0/PI)*m.log(R/r)
    # L = (1./2./pi).*log(R./r)
    return L
def serCoax_C(r,R):
    C = 2.0*PI*SI.EPS0/m.log(R/r)
#   C = 2.*pi./log(R./r)
    return C
def coax_Zw_LC(L,C):
    Zw=m.sqrt(L/C)
    return Zw # Ohm
def coax_Zw(r,R):
    Zw=coax_Zw_LC(serCoax_L(r,R), serCoax_C(r,R))
    return Zw # Ohm
def coax_Zw2(r,R,eps):
    Zw=(60.0/m.sqrt(eps))*m.log(R/r)
    return Zw # Ohm
def coaxGamma( f, r, R ):
    seriec_L = serCoax_L(r,R)
    series_C = serCoax_C(r,R)
    gamma = cm.sqrt( (0+1j*2*PI*f*seriec_L)*(0+1j*2*PI*f*series_C) )
    return gamma # Ohm
def coaxGamma_LC( f, seriec_L, series_C ):
    gamma = cm.sqrt( (0+1j*2*PI*f*seriec_L)*(0+1j*2*PI*f*series_C) )
    return gamma # Ohm
def coaxLambda( coaxGamma ):
    lambda_coax = 2*PI/imag(coaxGamma)
    return lambda_coax # m

def GammaTRL( f, seriec_L, series_C, seriec_R, series_G ):
    w = 2*PI*f
    gamma = cm.sqrt( ( seriec_R+1j*w*seriec_L )*( series_G + 1j*w*series_C ) )
    Zw = cm.sqrt( (seriec_R + 1j*w*seriec_L)/(series_G + 1j*w*series_C) )
    return gamma, Zw # m



