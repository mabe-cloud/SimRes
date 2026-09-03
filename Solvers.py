import numpy as np
import pandas as pd

#Class LinearSystemSolve():


def TDMA(T_matriz, D):
    a = np.diagonal(T_matriz, offset = - 1)
    b = np.diagonal(T_matriz, offset=0)
    c = np.diagonal(T_matriz, offset=+1)
    d = D

    n = len(d)
    c_ = np.zeros(n-1)
    d_ = np.zeros(n)
    x = np.zeros(n)

    # Forward elimination
    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]
    for i in range(1, n-1):
        c_[i] = c[i] / (b[i] - a[i-1] * c_[i-1])
    for i in range(1, n):
        d_[i] = (d[i] - a[i-1] * d_[i-1]) / (b[i] - a[i-1] * c_[i-1])

    # Back substitution
    x[-1] = d_[-1]
    for i in range(n-2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i+1]

    return x

# Critério de Scarvorought, 1966
n = 6 # Números de algarismos significativos
Eppara = 0.5*10**(2-n) # Termo relativos


def Gauss_Seidel(A, B, x0):
    x = np.array(x0)
    iter = 0
    n = []
    estimativa = []
    E_pest = []
    Epest = np.linspace(100, 100, len(B))
    v_iter = []

    while np.max(Epest) >= Eppara(6) and iter <= 1000:
        x_old = np.copy(x)

        for i in range(len(B)):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i + 1:], x_old[i + 1:])
            x[i] = (B[i] - sum1 - sum2) / A[i, i]

        Epest = np.abs((x - x_old) / x) * 100

        iter += 1

        v_iter.append(iter)
        estimativa.append(x.copy())
        E_pest.append(Epest)

    return x, v_iter, estimativa, E_pest


def GaussPiv(A, B):
    A = A.astype(float)
    B = B.astype(float)

    S = np.hstack((A, B.reshape(-1, 1)))

    n = A.shape[0]
    for k in range(n - 1):
        ipr = np.argmax(np.abs(S[k:, k])) + k

        if ipr != k:
            S[[k, ipr], :] = S[[ipr, k], :]

        for i in range(k + 1, n):
            fator_m = S[i, k] / S[k, k]
            S[i, k:] -= fator_m * S[k, k:]

    X = np.zeros(n)
    X[-1] = S[-1, -1] / S[-1, -2]
    for i in range(n - 2, -1, -1):
        X[i] = (S[i, -1] - np.dot(S[i, i + 1:n], X[i + 1:])) / S[i, i]

    return X


def Gauss_Seidel_rlx(A, B, l, x0):
    x = np.array(x0)
    iter = 0
    n = []
    estimativa = []
    E_pest = []
    Epest = np.linspace(100, 100, len(B))
    v_iter = []

    while np.max(Epest) >= Eppara(6) and iter <= 1000:
        x_old = np.copy(x)

        for i in range(len(B)):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i + 1:], x_old[i + 1:])
            x[i] = (B[i] - sum1 - sum2) / A[i, i]
            x[i] = l * x[i] + (1 - l) * x_old[i]

        # Critério de parada
        Epest = np.abs((x - x_old) / x) * 100

        iter += 1

        # Salvando o número de iterações e as estimativas dos erros
        v_iter.append(iter)
        estimativa.append(x.copy())
        E_pest.append(Epest)

    return x, v_iter, estimativa, E_pest