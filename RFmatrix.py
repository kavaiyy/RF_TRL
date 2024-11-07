import cmath as cm
import CircElements as CircElements
import numpy as np


def reverse_A_2x2(A):
    Ar = np.zeros( (2,2), dtype='complex' )
    A_det = A[0,0]*A[1,1]-A[0,1]*A[1,0]
    Ar[0,0]=A[1,1]/A_det
    Ar[1,1]=A[0,0]/A_det
    Ar[0,1]=-A[0,1]/A_det
    Ar[1,0]=-A[1,0]/A_det
    return 


def ABCD_TRL( Z, gamma, Lth ):
    A = np.zeros((2,2), dtype='complex')
    A[0,0] = cm.cosh(gamma*Lth)         # A
    A[0,1] = Z*cm.sinh(gamma*Lth)       # B
    A[1,0] = (1/Z)*cm.sinh(gamma*Lth)   # C
    A[1,1] = cm.cosh(gamma*Lth)         # D
    return A

def A_par_Z(Z):
    A = np.zeros((2,2), dtype='complex')
    A[0,0] = 1.0
    A[0,1] = 0.0
    A[1,0] = 1/Z
    A[1,1] = 1.0
    return A

def A_ser_Z(Z):
    A = np.zeros((2,2), dtype='complex')
    A[0,0] = 1.0
    A[0,1] = Z
    A[1,0] = 0.0
    A[1,1] = 1.0
    return A


def A_stub( Zw, gamma, Lth ):
    # A = np.zeros((2,2), dtype='complex')
    Z_stub = CircElements.Z_stub(Zw, gamma, Lth)
    return A_par_Z(Z_stub)
# #  Capacitance in the end of Stub
def A_stub_C_load(f,Zw,gamma, Lth, C):
    Z = CircElements.Z_stub_C_load(f,Zw,gamma, Lth, C);
    return A_par_Z(Z)
# #  Capacitance in the middle of Stub
def A_stub_C(f,Zw,gamma, Lth_1, C, Lth_2):
    Z = CircElements.Z_stub_C(f,Zw,gamma, Lth_1, C, Lth_2);
    return A_par_Z(Z)
def A_stub_Open(Zw, gamma, Lth):
    Z = CircElements.Z_stub_open(Zw,gamma,Lth);
    return A_par_Z(Z)


# # #
def A_P_equiv_scheme(Y1,Y2,Y3):
    A = np.zeros((2,2), dtype='complex')
    A[0,0] = 1.0 + Y2/Y3
    A[0,1] = 1/Y3
    A[1,0] = Y1 + Y2 + (Y1*Y2)/Y3
    A[1,1] = 1.0 + Y1/Y3
    return A
# # #
def A_T_equiv_scheme(Z1,Z2,Z3):
    A = np.zeros((2,2), dtype='complex')
    A[0,0] = 1.0 + Z1/Z3
    A[0,1] = Z1 + Z2 + (Z1*Z2)/Z3
    A[1,0] = 1/Z3
    A[1,1] = 1.0 + Z2/Z3
    return A


# # NO
def A_stub_exp( Zw, gamma, Lth ):
    Z = CircElements.Z_stub(Zw, gamma, Lth)
    return A_ser_Z(Z)
# # NO

# # # Zin for TRL line loaded with Z on the end

def Zin_ABCD(A, Z_load):
    Zin = ( A[0,0]*Z_load+A[0,1] )/( A[1,0]*Z_load+A[1,1] )
    return Zin
# Additional in case:
def G_refl_coef(Zin, Zw):
    return ( Zin - Zw )/( Zin + Zw )
def ABCD_load(A, Z_load,Zw):
    Zin = ( A[0,0]*Z_load+A[0,1] )/( A[1,0]*Z_load+A[1,1] )
    return Zin, G_refl_coef(Zin, Zw)

# # # S-parameters functions:
def Spar_passivity(S):
    return abs(S[0,0])**2 + abs(S[1,0])**2
def S_equation_resonance(S):
    return abs( S[1,1]-S[1,0]*(S[0,1]/S[0,0]) )



# # # By David M Pozar
def ABCD_to_Z(A):
    A,B,C,D = Am[0,0], Am[0,1], Am[1,0],Am[1,1]
    A_det = np.linalg.det(A)
    Z = np.zeros( (2,2), dtype='complex')
    Z[0,0]=A/C
    Z[0,1]=A_det/C
    Z[1,0]=1/C
    Z[1,1]=D/C
    return Z
# # # By David M Pozar
def ABCD_to_Ypar(Am):
    A,B,C,D = Am[0,0], Am[0,1], Am[1,0],Am[1,1]
    Y = np.zeros( (2,2), dtype='complex')
    Y[0,0] = D/B
    Y[0,1] = ( B*C - A*D )/B
    Y[1,0] = -1.0/B
    Y[1,1] = A/B
    return Y
# # # By David M Pozar
def Ypar_to_ABCD(Y):
    ABCD = np.zeros( (2,2), dtype='complex')
    det = np.linalg.det(Y)
    ABCD[0,0] = -Y[1,1]/Y[1,0]
    ABCD[0,1] = -1.0/Y[1,0]
    ABCD[1,0] = -det/Y[1,0]
    ABCD[1,1] = -Y[0,0]/Y[1,0]
    return ABCD

# # # By David M Pozar: # # in matlab name as ABCD_to_Spar2()
def ABCD_to_Spar(Am,z):
    A,B,C,D = Am[0,0], Am[0,1], Am[1,0],Am[1,1]
    #
    S = np.zeros( (2,2), dtype='complex')
    #
    summ = A+B/z+z*C+D
    S[1,1] = ( A+B/z-z*C-D )/summ
    S[1,2] = 2*( A*D - B*C )/summ
    S[2,1] = 2.0/summ
    S[2,2] = (-A+B/z-z*C+D)/summ
    return S

# # # By David M Pozar:
def Transfer_Parameters(Am,z):
    # A,B,C,D = Am[0,0], Am[0,1], Am[1,0],Am[1,1]
    # #
    # S = np.zeros( (2,2), dtype='complex')
    # #
    # summ = A+B/z+z*C+D
    # S[1,1] = ( A+B/z-z*C-D )/summ
    # S[1,2] = 2*( A*D - B*C )/summ
    # S[2,1] = 2.0/summ
    # S[2,2] = (-A+B/z-z*C+D)/summ
    return 0



def Smatrix_Hybrid_90degrees():
    Am = np.zeros( (4,4), dtype='complex')
    Am[0,1],Am[1,0], Am[2,3],Am[3,2] = 1.0j,1.0j, 1.0j,1.0j
    Am[0,2],Am[2,0], Am[1,3],Am[3,1] = 1.0,1.0, 1.0,1.0
    return Am*(-1/2**0.5)




