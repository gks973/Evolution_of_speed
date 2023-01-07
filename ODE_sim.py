import numpy as np
import random
import matplotlib.pyplot as plt

## Introducing all strains simultaneously

dt = 0.1
# Runtime of simulation:
T = int(round(500000/dt))
#T = int(round(10000/dt))

Nv = 6 # Number of variants

I0 = 0.001

Ts = np.array([10, 50, 100, 200, 300, 400]) # I-state durations
taus = 1*Ts

#Ts = np.array([10, 10, 10, 10]) # I-state durations
#taus = np.array([10, 100, 200, 300])



S = np.zeros(T)
E = np.zeros(T)
I = np.zeros(T)
Q = np.zeros(T)
R = np.zeros(T)
S[0] = 1-I0

Iv = np.zeros((T,Nv))
for v in range(Nv):
    Iv[0,v]=I0/Nv

Ev = np.zeros((T,Nv))
Qv = np.zeros((T,Nv))

beta = 1 # Transmission rate
omega = 0.001 # Waning rate
q = 0.1

lambd = 1 # Recovery rate

def get_RI_t(Rv,t):
    T = Rv.shape[0]
    Nv = Rv.shape[1]
    R = np.sum(Rv[t,:])
    return R

def make_mutation(Iv_t, to_v, epsilon=0.001):
    r = np.random.rand()
    u = np.argmax(Iv_t) # The dominant variant
    
    v = to_v

    if Iv_t[u] <= epsilon:
        epsilon = 0.1 * Iv_t[u]
    if Iv_t[v] >= 1-epsilon:
        print("Not making mutation, recipient already present in large numbers")
    else:
        Iv_t[v] += epsilon
        Iv_t[u] -= epsilon
        print("Moving", epsilon, f"from {(u)} to {(v)}.")
    return Iv_t
    

ts = [0]
added = [0]
for t in range(1,T):
    ts.append(dt*t)
    E[t-1] = get_RI_t(Ev,t-1)
    I[t-1] = get_RI_t(Iv,t-1)
    Q[t-1] = get_RI_t(Qv,t-1)
    
    dE_dt_v = np.zeros((Nv))
    dI_dt_v = np.zeros((Nv))
    dQ_dt_v = np.zeros((Nv))
    dS_dt = omega * R[t-1]
    dR_dt = - omega * R[t-1]
    for v in range(Nv):
        dS_dt += -beta * S[t-1] * Iv[t-1,v]
        dE_dt_v[v] = beta * S[t-1] * Iv[t-1,v] - (1/taus[v]) * Ev[t-1,v]
        dI_dt_v[v] = (1/taus[v]) * Ev[t-1,v] - (1/Ts[v]) * Iv[t-1,v] - q * Iv[t-1,v]
        dQ_dt_v[v] = q * Iv[t-1,v] - (1/Ts[v]) * Qv[t-1,v]
        dR_dt += (1/Ts[v]) * Iv[t-1,v] + (1/Ts[v]) * Qv[t-1,v]
    S[t] = S[t-1] + dt * dS_dt
    Ev[t,:] = Ev[t-1,:] + dt * dE_dt_v
    Iv[t,:] = Iv[t-1,:] + dt * dI_dt_v
    Qv[t,:] = Qv[t-1,:] + dt * dQ_dt_v
    R[t] = R[t-1] + dt * dR_dt
    
    if t % 10000 == 0:
        print(t, "out of", T)
        
#%matplotlib inline

#%%

from matplotlib.ticker import FuncFormatter

fontsize=20

font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : fontsize}

mystyle = 'default' #'seaborn'

plt.style.use(mystyle)
plt.rc('font', **font)


fig, axes= plt.subplots(1, 1, dpi=175, figsize=(8,6))

#plt.plot(ts, S, label="S", alpha=0.5)
#plt.plot(ts, I, label=r"$I = \sum_n I_n$", alpha=0.5)

for v in range(Nv):
    plt.plot(ts, Ev[:,v]+Iv[:,v]+Qv[:,v], "-", label=f"T={Ts[v]}")
#plt.plot(ts, R, label=r"$R = \sum_n R_n$", alpha=0.5)
#plt.plot(ts, S+I+R, label="S + I + R")
#plt.ylim([0,0.3])
plt.xlim([1,T*dt])
plt.legend(frameon=True,loc="lower right")
#plt.title(f"q={q}. All variants present initially.")
plt.xlabel("Time (Days)")
plt.ylabel("Total exposed and infected (E+I+Q)")
plt.xscale('log')

plt.show()
