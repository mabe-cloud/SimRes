import numpy as np
from scipy.special import erfc

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

