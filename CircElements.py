import cmath as cm

PI = cm.pi
def coth(x):
    return 1/cm.tanh(x)


def Z_C(f,C):
    w = 2*PI*f
    return -1j*1/(w*C)
def Z_L(f,L):
    w = 2*PI*f
    return w*L*1j
def Z_real_C(f,C, Rs,Li,Rd): # Rd should be very large
    w = 2*PI*f
    Z = 1j*w*Li + Rs + Rd*Z_C(f,C)/(Rd+Z_C(f,C))
    return Z



def Z_stub(Zw,gamma,length):
    Z = Zw*cm.tanh(gamma*length)
    return Z
def Y_stub(Zw,gamma,length):
    Z_input = CircElements.Z_stub(Zw,gamma,length)
    Y = 1/Z_input
    return Y



def Z_stub_open(Zw,gamma,length):
    Z = Zw*coth(gamma*length)
    return Z






