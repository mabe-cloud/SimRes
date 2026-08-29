import numpy as np
from scipy.special import erfc
from scipy.special import expi

# Escoamento Linear

def p1_1D(x, t, pe, pw, L, k, phi, mu, ct, N=100):
    eta = k / (phi * mu * ct)
    soma = np.zeros_like(x, dtype=float)
    
    for n in range(1, N + 1):
        a = (n * np.pi) / L
        termo = (np.exp(-(a**2)*eta*t)/n)*np.sin(a*x)
        soma += termo

    p = (pe - pw)*((x/L) + (2/np.pi)*soma) + pw
    
    return p

def p2_1D(x, t, p0, qw, mu, L, k, A, phi, ct):
    eta = k / (phi * mu * ct)
    a = (qw * mu * L)/(k*A)
    b = 4*eta*t

    p = p0 - a*(np.sqrt(b/(np.pi*L**2))*np.exp(-(x**2)/b) - ((x/L)*erfc(x/np.sqrt(b))))
    
    return p

def p_transiente_1D_radial(r, t, p0, qw, mu, h, k, phi, ct):
    a = (qw *mu)/(4 * np.pi * k * h)
    b = (phi * mu * ct * r**2)/(4 * k * t)
    return p0 + a*expi(-b)

def p_pseudopermanente_1D_radial(r, re, rw, t, p0, qw, mu, h, k, phi, ct):
    a = (qw *mu)/(2 * np.pi * k * h)
    b = 2*k*t/(phi * mu * ct * re**2)
    return p0 - a * (b - np.log(r/rw) + 1/2*(r/re)**2 + np.log(re/rw) - 3/4)

