import numpy as np
import matplotlib.pyplot as plt
from analiticas import p1_1D, p2_1D, p_pseudopermanente_1D_radial, p_transiente_1D_radial

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
    plt.plot(x, p, label=f"t = {t} s")

plt.xlabel('Posição (m)', size=13)
plt.ylabel("Pressão (Pa)", size=13)
plt.title("Solução Regime Permanente Linear", size=16)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

#%% Regime Transiente ---------------------------------------------------------

tempos = [1, 2, 3, 4, 5]    

plt.figure(figsize=(9, 6))

for t in tempos:
    p = p2_1D(x, t,p0, qw, mu, L, k, A, phi, ct)
    plt.plot(x, p, label=f"t = {t} s")

plt.xlabel('Posição (m)', size=13)
plt.ylabel("Pressão (Pa)", size=13)
plt.title("Regime Transiente Linear", size=16)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

#%% Pseudopermanente ----------------------------------------------------------

tempos_longos = [10, 15, 20, 25, 30, 40, 50, 70, 100]

plt.figure(figsize=(10, 7))

for t in tempos_longos:
    p = p2_1D(x, t, p0, qw, mu, L, k, A, phi, ct)
    plt.plot(x, p, label=f"t = {t} s")

plt.axvline(x=L,color="black",linestyle="--",alpha=0.5)


plt.xlabel('Posição (m)', size=13)
plt.ylabel("Pressão (Pa)", size=13)
plt.title("Regime Pseudopermanente Linear", size=16)
plt.grid(True, alpha=0.3)
plt.legend(ncol=2)
plt.tight_layout()
plt.show()


####################################
############## Radial ##############
####################################

# Dados

pe = 300 * 98066.5 # kgf/cm2 para Pa
pw = 150 * 98066.5 # kgf/cm2 para Pa
rw = 1 # m
re = 500 # m
k = 20 * 9.869e-16 # md para m2
phi = 0.18
mu = 0.8 *1e-3 # cp para Pa.s
ct = 150e-6 * 1.0197e-5 # (kgf/cm2)^-1 para Pa^-1
N = 100
h = 20 # espessura da formacao - m

# Dados para o transiente

Bo = 1.2 # m3 / m3std
qw_std = 400 * 1/86400 # m3std/dia para m3std/s
p0 = pe # pressão inicial
qw = Bo * qw_std

r = np.linspace(rw, re, 500)

plt.figure(figsize=(9, 6))
tempos = [60*30, 60*60*3, 60*60*12, 86400*2, 4*86400, 8*86400]
for t in tempos:
    p = p_transiente_1D_radial(r, t, p0, qw, mu, h, k, phi, ct) * 1e-6
    plt.plot(r, p, label=f"t = {t} s = {t/(60*60)} h")

plt.xlabel('Posição radial (m)', size=13)
plt.ylabel("Pressão (MPa)", size=13)
plt.title(r'Solução Transiente Radial 1D',size=16)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

tempos = [60*30, 86400, 86400*3, 7*86400, 14*86400, 30*86400]
for t in tempos:
    p = p_pseudopermanente_1D_radial(r, re, rw, t, p0, qw, mu, h, k, phi, ct) * 1e-6
    plt.plot(r, p, label=f"t = {t} s = {t/(60*60)} h")

plt.xlabel('Posição radial (m)', size=13)
plt.ylabel("Pressão (MPa)", size=13)
plt.title(r'Solução Pseudopermanente Radial 1D',size=16)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
