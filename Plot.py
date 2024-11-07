import numpy as np
import math as m

# import mypackage.someFuncs as myFunction
# import MeshGrids as MeshGrids


import CircElements as RF_elements
import CoaxLine as Coax_TRL
import ICRH_2D_Ant as Plasm_2D
import RFmatrix
cm = 10.0**-2

mF  = 10**-3
mkF = 10**-6
nF  = 10**-9
pF  = 10**-12

mH  = 10**-3
mkH = 10**-6
nH  = 10**-9
pH  = 10**-12


C = 10.0 * nF
L = 10.0 * mkH

# α = 31.0
# print(α)
# Δ  = 33.0
# print(Δ)


# # *********************************************************************
# # *********************** Grid definition ***********************
f_MHz = np.arange(1.0,300.0)#MHz
f = f_MHz*10**6
Z_C = np.zeros((len(f_MHz)), dtype='complex')
Z_L = np.zeros((len(f_MHz)), dtype='complex')
Z_stub = np.zeros((len(f_MHz)), dtype='complex')
Y_stub = np.zeros((len(f_MHz)), dtype='complex')
Y_stub_C = np.zeros((len(f_MHz)), dtype='complex')
Y_stub_midC = np.zeros((len(f_MHz)), dtype='complex')
Y_stub_Z = np.zeros((len(f_MHz)), dtype='complex')
Y_stub_midZ = np.zeros((len(f_MHz)), dtype='complex')
Z_ICRH_ant = np.zeros((len(f_MHz)), dtype='complex')
ABCD_plasm =np.zeros((len(f_MHz),2,2), dtype='complex')
# # *********************** Grid definition ***********************
# # *********************************************************************
# # *********************************************************************
# # *********************** Y axis definition ***********************


# # *********************** Y axis definition ***********************
# # *********************************************************************
# Z_stub(Zw,gamma,length)
# Y_stub
Len1  = 1.0 	 # m


# # ------------------------------------------------
# # Coaxial Line Parameters:
r_in  = 2.75 *cm # m
R_out = 6.0  *cm # m
Zw = Coax_TRL.coax_Zw( r_in, R_out )
# # ------------------------------------------------


# # ------------------------------------------------
# # plasm Parameters:
# L_plasm = 35.0*cm # 35*cm %19*cm
# N0 = 18.0 * 10.0**18 # % m^-3 L-mode: N0 = 1.0 * 10.^18; dx = 4.8*cm;
# dx = 8.8*cm  #1*cm
# w = 11.8*cm  #
# wz = w/2     #
# s = 1.6*cm   #
# d = 3.8*cm   #
# t = 0.94*cm  #
ant_N0,ant_dx,ant_w,ant_wz,ant_s,ant_d,ant_t,ant_Lplasm = Plasm_2D.def_params()
ant_Lplasm = 35.0*cm
w_ci = Plasm_2D.plasmParam_GlobusM1()
serPlasm_C = Plasm_2D.AntennaICRH__ser_C(ant_s,ant_d,ant_wz,ant_t)
# # ------------------------------------------------
def get_plasm( f ):
    alph = Plasm_2D.alpha( f, ant_N0,  ant_dx,  w_ci )
    N = Plasm_2D.N_(alph, ant_dx, ant_d)
    del11 = Plasm_2D.delta(alph,ant_dx,ant_d)
    serPlasm_L = Plasm_2D.AntennaICRH__ser_L( ant_wz,ant_d,N,del11 )
    serPlasm_R = Plasm_2D.AntennaICRH__ser_R( f, ant_wz, ant_d, alph, del11 )
    gamma_plasm, plasm_Zw = Coax_TRL.GammaTRL( f, serPlasm_L, serPlasm_C, serPlasm_R, 0 );
    return gamma_plasm, plasm_Zw

    # % *************** the plasma ******
    # alph = ICRH_2D_Ant.alpha( f,N0,dx,w_ci );
    # N = ICRH_2D_Ant.N_(alph,dx,d);
    # del = ICRH_2D_Ant.delta(alph,dx,d);
    # serPlasm_L = ICRH_2D_Ant.series_L( wz,d,N,del );
    # serPlasm_R = ICRH_2D_Ant.series_R( f, wz, d, alph, del );
    # plasm_ser_L = serPlasm_L;
    # plasm_ser_R = serPlasm_R;
    # [gamma, plasm_Zw] = CoaxLine.GammaTRL( f, serPlasm_L, serPlasm_C, serPlasm_R, 0 );
    # gamma_Plasm = gamma;
    # %      arr_plasm_Zw(j) = plasm_Zw;
    # Aplsm(:,:,j) = RFmatrix.ABCD_TRL( plasm_Zw, gamma, L_plasm );
    # % *************** the plasma ******


print("start cycle:")
for i in range(0, len(f_MHz)):
	gamma_plasm, Zw_plasm = get_plasm( f[i] )
	ABCD_plasm[i,:,:] = RFmatrix.ABCD_TRL( Zw_plasm, gamma_plasm, ant_Lplasm )
	Z_ICRH_ant[i] = RFmatrix.Zin_ABCD(ABCD_plasm[i,:,:], 0.0j)
	# Z_C[i] = RF_elements.Z_C(f[i], C)
	# Z_L[i] = RF_elements.Z_L(f[i], L)
	gamma = Coax_TRL.coaxGamma( f[i], r_in, R_out )
	Z_stub[i] = RF_elements.Z_stub(Zw, gamma, Len1)
	# Y_stub[i] = RF_elements.Y_stub(f[i], L)
	# Y_stub_C[i] = RF_elements.Z_L(f[i], L)
	# Y_stub_midC[i] = RF_elements.Z_L(f[i], L)
	# Y_stub_Z[i] = RF_elements.Z_L(f[i], L)
	# Y_stub_midZ[i] = RF_elements.Z_L(f[i], L)
	pass
print("end cycle:")






import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mypackage.CGS import inch_to_cm
# plt.style.use("MyClassic.mplstyle")
plt.style.use("MyClassic_From_stylelib") # "ggplot" "MyClassic" "classic" "MyClassic_From_stylelib"

fig = plt.figure(figsize=(8.0, 6.0))

# plt.plot(f_MHz, Z_C.imag, color="crimson", linewidth=1.5, label="Z_C")
# plt.plot(f_MHz, Z_L.imag, color="blue", linewidth=1.5, label="Z_L")
# plt.plot(f_MHz, Z_stub.imag, color="blue", linewidth=1.5, label="Z_stub")
# plt.plot(f_MHz, 1/Z_stub.imag, color="blue", linewidth=1.5, label="Y_stub")


plt.plot(f_MHz, Z_ICRH_ant.real, color="blue", linewidth=2.5, label="ICRH_ant_Z_real")
plt.plot(f_MHz, Z_ICRH_ant.imag, color="red",  linewidth=1.0, label="ICRH_ant_Z_imag")




plt.legend(loc="upper right")
plt.xlabel("f, MHz")
plt.ylabel("Z, Ohm")






# plt.savefig("SomePlot.png",format="png", transparent=False, bbox_inches='tight')
plt.show()
# np.savetxt("meshGrid_mNz.txt", mNz)
# np.savetxt("meshGrid_mNy.txt", mNy)







# px = 1/plt.rcParams['figure.dpi']  # pixel in inches

# plt.yticks(fontsize=18)
# plt.xticks(fontsize=18)
# plt.xlabel("Ny", fontsize=20)
# plt.ylabel("Sp", fontsize=20)


# plt.axvline(x = nz1, color = 'black')
# plt.axvline(x = nz2, color = 'black')

# plt.axvline(x = ny1, color = 'black')
# plt.axvline(x = ny2, color = 'black')