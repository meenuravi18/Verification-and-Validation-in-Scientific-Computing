import numpy as np
import wildfire
from utils.functions import G
import matplotlib.pyplot as plt

kap = 0.1
eps = 0.3
upc = 3
alp = 0
q = 1
x_min, x_max = 0, 1
y_min, y_max = 0, 1
t_min, t_max = 0, 1

physical_parameters = {    
    'kap': kap, # diffusion coefficient
    'eps': eps, # inverse of activation energy
    'upc': upc, # u phase change
    'q': q, # reaction heat
    'alp': alp, 
    'x_lim': (x_min, x_max), # x-axis domain 
    'y_lim': (y_min, y_max), # y-axis domain
    't_lim': (t_min, t_max), # time domain
    'components':(True, False, True) #diffusion, convection, reaction
}


wildfire_ = wildfire.Fire(**physical_parameters)
h=[]
fuels=[]
temps=[]
Nx_vals = np.array([16, 32, 64, 128])
for Nx in Nx_vals:
    
    Ny = Nx
    Nt = 5 * Nx**2
    u0=lambda X, Y: np.exp(-((X-0.5)**2 + (Y-0.5)**2)/0.02) 
    b0 = lambda X, Y: np.ones_like(X)
    V=lambda x, y, t: (0*x, 0*y)
    t, X, Y, U, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)

    SRQ_betamax_num   = np.max(B)
    fuels.append(SRQ_betamax_num)
    SRQ_umax_num   = np.max(U)
    temps.append(SRQ_umax_num)
    h.append((x_max-x_min)/(Nx-1))
   
    
    print("done with Nx =", Nx)

print(fuels)
print(temps)
print("*****")
temp_triplet1cm=temps[0]-temps[1]
temp_triplet1mf=temps[1]-temps[2]
temp_triplet2cm=temps[1]-temps[2]
temp_triplet2mf=temps[2]-temps[3]
p1=np.log(temp_triplet1cm/temp_triplet1mf)/np.log(2)
p2=np.log(temp_triplet2cm/temp_triplet2mf)/np.log(2)
print(p1)
print(p2)
fuel_triplet1cm=fuels[0]-fuels[1]
fuel_triplet1mf=fuels[1]-fuels[2]
fuel_triplet2cm=fuels[1]-fuels[2]
fuel_triplet2mf=fuels[2]-fuels[3]
p3=np.log(fuel_triplet1cm/fuel_triplet1mf)/np.log(2)
p4=np.log(fuel_triplet2cm/fuel_triplet2mf)/np.log(2)
print(p3)
print(p4)
print("*****")
#factor of safety check
print("*******************Medium to Fine*******************")
fs=0
if abs((p2-2)/2) <= 0.1:
    fs=1.25
else:
    fs=3.0

print("factor of safety:", fs)
gci_temp=((fs)/((2**p2)-1))*(abs(temp_triplet2mf))
print("GCI for temperature:", gci_temp)
# fuel
fs=0
if abs((p4-2)/2) <= 0.1:
    fs=1.25
else:
    fs=3.0
print("factor of safety:", fs)
gci_fuel=((fs)/((2**p4)-1))*(abs(fuel_triplet2mf))
print("GCI for fuel:", gci_fuel)



de=(temps[2]-temps[3])/((2**p2)-1)
print("de for temperature:", de)

de=(fuels[2]-fuels[3])/((2**p4)-1)
print("de for fuel:", de)



print("*******************COARSE TO MEDIUM*******************")
fs=0
if abs((p1-2)/2) <= 0.1:
    fs=1.25
else:
    fs=3.0

print("factor of safety:", fs)
gci_temp=((fs)/((2**p1)-1))*(abs(temp_triplet1cm))
print("GCI for temperature:", gci_temp)
# fuel
fs=0
if abs((p3-2)/2) <= 0.1:
    fs=1.25
else:
    fs=3.0
print("factor of safety:", fs)
gci_fuel=((fs)/((2**p3)-1))*(abs(fuel_triplet1cm))
print("GCI for fuel:", gci_fuel)

ubar = temps[2] + ((temps[2] - temps[1])/(2**p1 - 1))
de = temps[1] - ubar
print("de for temperature:", de)

ubar = fuels[2] + ((fuels[2] - fuels[1])/(2**p3 - 1))
de = fuels[1] - ubar
print("de for fuel:", de)






# Work for Round off error calculation
# Run simulation with desired precision
# Re-run simulation with higher precision
# Estimate error in desired precision simulation
# use same mesh/time step for each
# Iteratively converge each case down to machine precision (i.e., until iterative
# residuals can no longer be reduced due to round-off error)

# I use 16 as the mesh size
# size=128
# kap = np.float32(0.1)
# eps = np.float32(0.3)
# upc = np.float32(3)
# alp = np.float32(0)
# q = np.float32(1)
# x_min, x_max = np.float32(0), np.float32(1)
# y_min, y_max = np.float32(0), np.float32(1)
# t_min, t_max = np.float32(0), np.float32(1)

# physical_parameters = {    
#     'kap': kap, # diffusion coefficient
#     'eps': eps, # inverse of activation energy
#     'upc': upc, # u phase change
#     'q': q, # reaction heat
#     'alp': alp, 
#     'x_lim': (x_min, x_max), # x-axis domain 
#     'y_lim': (y_min, y_max), # y-axis domain
#     't_lim': (t_min, t_max), # time domain
#     'components':(True, False, True) #diffusion, convection, reaction
# }


# wildfire_ = wildfire.Fire(**physical_parameters)
# h=[]
# fuels=np.float32(0)
# temps=np.float32(0)
# print(temps.dtype)
# Nx_vals = np.array([size])
# for Nx in Nx_vals:
    
#     Ny = Nx
#     Nt = 5 * Nx**2
#     u0=lambda X, Y: np.float32(np.exp(-((X-0.5)**2 + (Y-0.5)**2)/0.02)) 
#     b0 = lambda X, Y: np.float32(np.ones_like(X))
#     V=lambda x, y, t: (0*x, 0*y)
#     t, X, Y, U, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)

#     SRQ_betamax_num   = np.max(B)
#     fuels=np.float32(SRQ_betamax_num)
#     SRQ_umax_num   = np.max(U)
#     temps=np.float32(SRQ_umax_num)
#     print("done with Nx =", Nx)
# #####################################################################################################################

# kap = np.longdouble(0.1)
# eps = np.longdouble(0.3)
# upc = np.longdouble(3)
# alp = np.longdouble(0)
# q = np.longdouble(1)
# x_min, x_max = np.longdouble(0), np.longdouble(1)
# y_min, y_max = np.longdouble(0), np.longdouble(1)
# t_min, t_max = np.longdouble(0), np.longdouble(1)

# physical_parameters = {    
#     'kap': kap, # diffusion coefficient
#     'eps': eps, # inverse of activation energy
#     'upc': upc, # u phase change
#     'q': q, # reaction heat
#     'alp': alp, 
#     'x_lim': (x_min, x_max), # x-axis domain 
#     'y_lim': (y_min, y_max), # y-axis domain
#     't_lim': (t_min, t_max), # time domain
#     'components':(True, False, True) #diffusion, convection, reaction
# }
# fuels_high=np.longdouble(0)
# temps_high=np.longdouble(0)
# print(temps_high.dtype)
# wildfire_ = wildfire.Fire(**physical_parameters)
# Nx_vals = np.array([size])
# for Nx in Nx_vals:
    
#     Ny = Nx
#     Nt = 5 * Nx**2
#     u0=lambda X, Y: np.longdouble(np.exp(-((X-0.5)**2 + (Y-0.5)**2)/0.02))
#     b0 = lambda X, Y: np.longdouble(np.ones_like(X))
#     V=lambda x, y, t: (0*x, 0*y)
#     t, X, Y, U, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)

#     SRQ_betamax_num   = np.max(B)
#     fuels_high=np.longdouble(SRQ_betamax_num)
#     SRQ_umax_num   = np.max(U)
#     temps_high=np.longdouble(SRQ_umax_num)
#     print("done with Nx =", Nx)

# roundoff_fuel = (np.abs(fuels - fuels_high))
# roundoff_temp = (np.abs(temps - temps_high))
# print(temps.dtype)
# print(temps_high.dtype)

# print("Round-off fuel:", roundoff_fuel)
# print("Round-off temperature:", roundoff_temp)
