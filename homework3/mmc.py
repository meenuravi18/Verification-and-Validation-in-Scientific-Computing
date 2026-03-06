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
    'components':(True, False, False)
}


wildfire_ = wildfire.Fire(**physical_parameters)
u_exacts=[]
h=[]
srq_errors=[]
l2normerrors = []
linfnormerrors=[]
Nx_vals = np.array([16, 32, 64, 128])
for Nx in Nx_vals:
    
    Ny = Nx
    Nt = 5 * Nx**2
    V = lambda x, y, t: (0*x, 0*y)
    b0 = lambda X, Y: np.sin(np.pi*X)*np.sin(np.pi*Y)
    u0 = lambda X, Y: np.sin(np.pi*X)*np.sin(np.pi*Y)

    t, X, Y, U, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)
    u_exact = (1 + t[-1]) * np.sin(np.pi*X) * np.sin(np.pi*Y)
    l2normerror = np.sqrt(((1/(Nx-1))**2) * (np.sum(np.abs(U - u_exact)**2)))
    linfnormerror = np.max(np.abs(U - u_exact))
    u_exacts.append(u_exact)
    l2normerrors.append(l2normerror) 
    linfnormerrors.append(linfnormerror)  

    SRQ_umax_num   = np.max(U)
    SRQ_umax_exact = np.max(u_exact)
    SRQ_umax_err   = abs(SRQ_umax_num - SRQ_umax_exact) 
    srq_errors.append(SRQ_umax_err)
    h.append((x_max-x_min)/(Nx-1))
    if Nx==128:
        error = np.abs(U - u_exact)
        plt.figure()
        plt.contourf(X, Y, error, 20)
        plt.colorbar(label="Error")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Local discretization error (Mesh Node=128) in space")
        plt.show()
    
    print("done with Nx =", Nx)

print(f"l2normerrors: {l2normerrors}")
print(f"linfnormerrors: {linfnormerrors}")
print(f"srq_errors: {srq_errors}")
h_arr = np.array(h)
E2 = np.array(l2normerrors)
Einf = np.array(linfnormerrors)
ESRQ = np.array(srq_errors)

p2 = np.log(E2[:-1]/E2[1:]) / np.log(h_arr[:-1]/h_arr[1:])
pinf = np.log(Einf[:-1]/Einf[1:]) / np.log(h_arr[:-1]/h_arr[1:])
psrq = np.log(ESRQ[:-1]/ESRQ[1:]) / np.log(h_arr[:-1]/h_arr[1:])

plt.figure()
plt.semilogx(h_arr[:-1], p2, 'o-', label='Order (L2)')
plt.semilogx(h_arr[:-1], pinf, 's-', label='Order (Linf)')
plt.semilogx(h_arr[:-1], psrq, '^-', label='Order (SRQ)')

plt.axhline(2, color='k', linestyle='--', label='2nd order')

plt.xlabel('h')
plt.ylabel('Observed Order p')
plt.legend()
plt.show()


# Work for Round off error calculation
# Run simulation with desired precision
# Re-run simulation with higher precision
# Estimate error in desired precision simulation
# use same mesh/time step for each
# Iteratively converge each case down to machine precision (i.e., until iterative
# residuals can no longer be reduced due to round-off error)

# I use 16 as the mesh size

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
#     'components':(True, False, False)
# }


# wildfire_ = wildfire.Fire(**physical_parameters)
# u_exacts=[]
# h=[]
# srq_errors=[]
# l2normerrors = []
# linfnormerrors=[]

# Nx_vals = np.array([16])
# for Nx in Nx_vals:
    
#     Ny = Nx
#     Nt = 5 * Nx**2 
#     V = lambda x, y, t: (0*x, 0*y)
#     b0 = lambda X, Y: (np.sin(np.pi*X)*np.sin(np.pi*Y)).astype(np.float32)
#     u0 = lambda X, Y: (np.sin(np.pi*X)*np.sin(np.pi*Y)).astype(np.float32)

#     t, X, Y, U_low, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)
#     U_low = U_low.astype(np.float32)
#     u_exact = (1 + t[-1]) * np.sin(np.pi*X) * np.sin(np.pi*Y)
#     l2normerror = np.sqrt(((1/(Nx-1))**2) * (np.sum(np.abs(U_low - u_exact)**2)))
#     linfnormerror = np.max(np.abs(U_low - u_exact))
#     SRQ_umax_num   = np.max(U_low)
#     SRQ_umax_exact = np.max(u_exact)
#     SRQ_umax_err   = abs(SRQ_umax_num - SRQ_umax_exact) 
   
#     h.append((x_max-x_min)/(Nx-1))
#     print("done with Nx =", Nx)


# # use long double precision for higher precision simulation
# kap = np.longdouble(0.1)
# eps = np.longdouble(0.3)
# upc = np.longdouble(3)
# alp = np.longdouble(0)
# q = np.longdouble(1)
# x_min, x_max = np.longdouble(0), np.longdouble(1)
# y_min, y_max = np.longdouble(0), np.longdouble(1)
# t_min, t_max = np.longdouble(0), np.longdouble(1)
# print(alp.dtype)

# physical_parameters = {    
#     'kap': kap, # diffusion coefficient
#     'eps': eps, # inverse of activation energy
#     'upc': upc, # u phase change
#     'q': q, # reaction heat
#     'alp': alp, 
#     'x_lim': (x_min, x_max), # x-axis domain 
#     'y_lim': (y_min, y_max), # y-axis domain
#     't_lim': (t_min, t_max), # time domain
#     'components':(True, False, False)
# }


# wildfire_ = wildfire.Fire(**physical_parameters)
# h=[]
# for Nx in Nx_vals:
    
#     Ny = Nx
#     Nt = 5 * Nx**2 
#     V = lambda x, y, t: (0*x, 0*y)
#     b0 = lambda X, Y: np.longdouble(np.sin(np.pi*X)*np.sin(np.pi*Y))
#     u0 = lambda X, Y: np.longdouble(np.sin(np.pi*X)*np.sin(np.pi*Y))

#     t, X, Y, U_high, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)

#     u_exact = (1 + t[-1]) * np.sin(np.pi*X) * np.sin(np.pi*Y)
#     l2normerror = np.sqrt(((1/(Nx-1))**2) * (np.sum(np.abs(U_high - u_exact)**2)))
#     linfnormerror = np.max(np.abs(U_high - u_exact))

#     SRQ_umax_num   = np.max(U_high)
#     SRQ_umax_exact = np.max(u_exact)
#     SRQ_umax_err   = abs(SRQ_umax_num - SRQ_umax_exact) 
   
#     h.append((x_max-x_min)/(Nx-1))
#     print("done with Nx =", Nx)


# roundoff_Linf = np.max(np.abs(U_low - U_high))
# roundoff_L2 = np.sqrt(((1/(Nx-1))**2) * np.sum((U_low - U_high)**2))
# print(U_low.dtype)
# print(U_high.dtype)

# print("Round-off Linf:", roundoff_Linf)
# print("Round-off L2:", roundoff_L2)
