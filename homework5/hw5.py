import numpy as np
import wildfire
from utils.functions import G
import matplotlib.pyplot as plt
from scipy.stats import qmc, norm
import random 

kap=0.1
eps = 0.3
upc = 4
alp = 0
q = 1
x_min, x_max = 0, 1
y_min, y_max = 0, 1
t_min, t_max = 0, 1
mu=3


def runSolver(physical_parameters):
    wildfire_ = wildfire.Fire(**physical_parameters)
    Nx=64   
    Ny = Nx
    Nt = 5 * Nx**2
    u0=lambda X, Y: np.exp(-((X-0.5)**2 + (Y-0.5)**2)/0.02) 
    b0 = lambda X, Y: np.ones_like(X)
    V=lambda x, y, t: (0*x, 0*y)
    t, X, Y, U, B = wildfire_.solvePDE(Nx, Ny, Nt, u0, b0, V, 'FD', 'RK4', last=True, acc=2, sparse=False)

    SRQ_betamax_num   = np.max(B)
    SRQ_umax_num   = np.max(U)

    print(f" SRQ_betamax_num is: {SRQ_betamax_num}")
    print(f" SRQ_umax_num is: {SRQ_umax_num}")

    print("done")
    return SRQ_umax_num

# Choosing an appropriate sigma value
# sigmas = [1.2, 1.3, 1.4]
# values=[]
# for i in sigmas:
#     print(f"Running for Sigma = {i}")
#     physical_parameters1 = {    
#         'kap': kap, # diffusion coefficient
#         'eps': eps, # inverse of activation energy
#         'upc': mu-i, # u phase change
#         'q':q, # reaction heat
#         'alp': alp, 
#         'x_lim': (x_min, x_max), # x-axis domain 
#         'y_lim': (y_min, y_max), # y-axis domain
#         't_lim': (t_min, t_max), # time domain
#         'components':(True, False, True) #diffusion, convection, reaction
#     }
#     physical_parameters2 = {    
#         'kap': kap, # diffusion coefficient
#         'eps': eps, # inverse of activation energy
#         'upc': mu+i, # u phase change
#         'q': q, # reaction heat
#         'alp': alp, 
#         'x_lim': (x_min, x_max), # x-axis domain 
#         'y_lim': (y_min, y_max), # y-axis domain
#         't_lim': (t_min, t_max), # time domain
#         'components':(True, False, True) #diffusion, convection, reaction
#     }
#     srq_high_u, srq_high_b = runSolver(physical_parameters2)
#     srq_low_u, srq_low_b = runSolver(physical_parameters1)

#     print("U difference:", abs(srq_high_u - srq_low_u))
#     print("B difference:", abs(srq_high_b - srq_low_b))


sigma=1.3
upcmin=mu-sigma
upcmax=mu+sigma
chi1=[0.55, 0.95, 1.0, 1.1, 1.5]
chi2=[0.1,0.4,0.6,0.75,0.8,0.9,0.91,0.97,1.3,1.6]
values=[]
for i in [upcmin, upcmax]:
    physical_parameters = {    
        'kap': kap, # diffusion coefficient
        'eps': eps, # inverse of activation energy
        'upc': i, # u phase change
        'q': q, # reaction heat
        'alp': alp, 
        'x_lim': (x_min, x_max), # x-axis domain 
        'y_lim': (y_min, y_max), # y-axis domain
        't_lim': (t_min, t_max), # time domain
        'components':(True, False, True) #diffusion, convection, reaction
    }
    values.append(runSolver(physical_parameters))
alpha=min(values)
beta=max(values)
expdata1=[alpha + (c*(beta-alpha)) for c in chi1]
expdata2=[alpha + (c*(beta-alpha)) for c in chi2]
print(expdata1)
print(expdata2)
#Part 2 Validation metric computation

# Initialize the LHS sampler (d = dimensions)

mean = mu        
std = sigma
l_bounds=0
# Transform to Normal Distribution (mean=0, std=1)
sampler10 = qmc.LatinHypercube(d=1, seed=42)
sample10= sampler10.random(n=10)
sample10 = sample10.flatten()
sample10 = norm.ppf(sample10, loc=mean, scale=std)

sampler25 = qmc.LatinHypercube(d=1, seed=43)
sample25= sampler25.random(n=25)
sample25 = sample25.flatten()
sample25 = norm.ppf(sample25, loc=mean, scale=std)

sampler100 = qmc.LatinHypercube(d=1, seed=44)
sample100= sampler100.random(n=100)
sample100 = sample100.flatten()
sample100 = norm.ppf(sample100, loc=mean, scale=std)
sample100[sample100 <= 0] = norm.ppf(np.random.random(np.sum(sample100 <= 0)), loc=mean, scale=std)

print(len(sample10), len(sample25), len(sample100))
print(np.min(sample10), np.min(sample25), np.min(sample100))

plt.figure()
plt.hist(sample10, bins=5, edgecolor='black')
plt.title("LHS Normal Samples, n=10")
plt.xlabel("upc")
plt.ylabel("Frequency")
plt.show()

plt.figure()
plt.hist(sample25, bins=8, edgecolor='black')
plt.title("LHS Normal Samples, n=25")
plt.xlabel("upc")
plt.ylabel("Frequency")
plt.show()

plt.figure()
plt.hist(sample100, bins=15, edgecolor='black')
plt.title("LHS Normal Samples, n=100")
plt.xlabel("upc")
plt.ylabel("Frequency")
plt.show()