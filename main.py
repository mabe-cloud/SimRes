import numpy as np
import matplotlib.pyplot as plt
from analiticas import p1_1D, p2_1D

p0 = 50000
pe = p0
pw = 1000
qw = 80
L = 1000
e = 0.2
A = L*e 
k = 0.01
phi = 0.20
mu = 0.4
ct = 1e-5

x = np.linspace(0, L, 500)


#%% Regime Permanente ---------------------------------------------------------
tempos = [50, 100, 200]
N = 100

plt.figure(figsize=(9, 6))

for t in tempos:
    p = p1_1D(x, t, pe, pw, L, k, phi, mu, ct, N)
    plt.plot(x, p, label=f"t = {t}")

plt.xlabel("x")
plt.ylabel("p")
plt.title("Solução Regime Permanente Linear")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

#%% Regime Transiente ---------------------------------------------------------

tempos = [1, 2, 3, 4, 5]    

plt.figure(figsize=(9, 6))

for t in tempos:
    p = p2_1D(x, t,p0, qw, mu, L, k, A, phi, ct)
    plt.plot(x, p, label=f"t = {t}")

plt.xlabel("x")
plt.ylabel("p")
plt.title("Regime Transiente Linear")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

#%% Pseudopermanente ----------------------------------------------------------

tempos_longos = [10, 15, 20, 25, 30, 40, 50, 70, 100]

plt.figure(figsize=(10, 7))

for t in tempos_longos:
    p = p2_1D(x, t, p0, qw, mu, L, k, A, phi, ct)
    plt.plot(x, p, label=f"t = {t}")

plt.axvline(x=L,color="black",linestyle="--",alpha=0.5)


plt.xlabel("x")
plt.ylabel("p")
plt.title("Regime Pseudopermanente Linear")
plt.grid(True, alpha=0.3)
plt.legend(ncol=2)
plt.tight_layout()
plt.show()
