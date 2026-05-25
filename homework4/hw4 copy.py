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


h=[]

Nx_vals = np.array([16, 32, 64, 128])
corner_points_kap=[0.08, 0.12]
corner_points_upc=[0.75, 1.75]
results = {}  
h = [(x_max - x_min) / (Nx - 1) for Nx in Nx_vals] 
for k in corner_points_kap:
    for u in corner_points_upc:
        fuels=[]
        temps=[]
        for Nx in Nx_vals:

            physical_parameters = {    
                'kap': k, # diffusion coefficient
                'eps': eps, # inverse of activation energy
                'upc': u, # u phase change
                'q': q, # reaction heat
                'alp': alp, 
                'x_lim': (x_min, x_max), # x-axis domain 
                'y_lim': (y_min, y_max), # y-axis domain
                't_lim': (t_min, t_max), # time domain
                'components':(True, False, True) #diffusion, convection, reaction
            }
            wildfire_ = wildfire.Fire(**physical_parameters)
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
        
            
            print("done with Nx =", Nx)
        results[(k, u)] = {'temps': temps, 'fuels': fuels}
for corner, data in results.items():
    k, u = corner
    print(f"\nWe are using the corner: kap={k}, upc={u}....................................")
    temps = data['temps']
    fuels = data['fuels']
    print("*****")
    print(fuels)
    print(temps)
    print("*****")
    temp_triplet1cm=temps[0]-temps[1]
    temp_triplet1mf=temps[1]-temps[2]
    temp_triplet2cm=temps[1]-temps[2]
    temp_triplet2mf=temps[2]-temps[3]
    p1=np.log(temp_triplet1cm/temp_triplet1mf)/np.log(2)
    p2=np.log(temp_triplet2cm/temp_triplet2mf)/np.log(2)
    print("====")
    print(p1)
    print(p2)
    print("====")
    fuel_triplet1cm=fuels[0]-fuels[1]
    fuel_triplet1mf=fuels[1]-fuels[2]
    fuel_triplet2cm=fuels[1]-fuels[2]
    fuel_triplet2mf=fuels[2]-fuels[3]
    p3=np.log(fuel_triplet1cm/fuel_triplet1mf)/np.log(2)
    p4=np.log(fuel_triplet2cm/fuel_triplet2mf)/np.log(2)
    print("*****")
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

all_gci_temp = []
all_gci_fuel = []

for corner, data in results.items():
    temps = data['temps']
    fuels = data['fuels']
    p2 = np.log((temps[1]-temps[2])/(temps[2]-temps[3])) / np.log(2)
    p4 = np.log((fuels[1]-fuels[2])/(fuels[2]-fuels[3])) / np.log(2)
    fs_t = 1.25 if abs((p2-2)/2) <= 0.1 else 3.0
    fs_f = 1.25 if abs((p4-2)/2) <= 0.1 else 3.0
    all_gci_temp.append((fs_t / (2**p2 - 1)) * abs(temps[2]-temps[3]))
    all_gci_fuel.append((fs_f / (2**p4 - 1)) * abs(fuels[2]-fuels[3]))

print(f"\nMax GCI temp across all corners: {max(all_gci_temp):.6f}")
print(f"Max GCI fuel across all corners: {max(all_gci_fuel):.6f}")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for corner, data in results.items():
    k, u = corner
    label = f'kap={k}, upc={u}'
    axes[0].plot(h, data['temps'], 'o-', label=label)
    axes[1].plot(h, data['fuels'], 's-', label=label)

axes[0].set_xlabel('h')
axes[0].set_ylabel('Max Temperature')
axes[0].set_title('Max Temperature vs h')
axes[0].legend()

axes[1].set_xlabel('h')
axes[1].set_ylabel('Max Fuel')
axes[1].set_title('Max Fuel vs h')
axes[1].legend()

plt.tight_layout()
plt.show()
