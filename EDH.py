import numpy as np
from scipy.special import erfc
from scipy.special import expi
from Solvers import Solve
import matplotlib.pyplot as plt
from scipy.stats import qmc, lognorm
import pandas as pd
import os

class Analytical:
    """
    Classe contendo as soluções analíticas.
    Suporta somente problemas 1D

    dimension : int
        Define a quantidade de dimensões espaciais do modelo.

    coordinates : 'Linear' ou 'Radial'
        Configura o tipo de modelagem a ser utilizada.

    ci : Pressão inicial
        Define a condição inicial do problema. É, estritamente, a pressão inicial no domínio do problema.


    cc : [[tipo da CC, tipo da CC],[lado esquerdo, lado direito]]
        São as condições de contorno do problema. Podem ser do tipo Dirichlet (pressão) ou Neumann (vazão).
        Considera-se que o poço está localizado no lado esquerdo e a fronteira do reservatório no lado direito.

    grid : [nx,ny,nz]
        Define o tamanho da malha do problema.

    system_units : 'SI', 'BR', 'USA' ou 'D'
        Informa o sistema de unidades a ser utilizado na modelagem, em que 'D' se refere a sem dimensão.

    """

    def __init__(self, dimension : int, coordinates : str, ci : float, cc : [list, list], grid : list, system_units = 'SI'):
        self.dimension, self.ci, self.cc, self.grid, self.system_units = dimension, ci, cc, grid, system_units
        self.coordinates = coordinates
        if any([None == var for var in self.__dict__.values()]):
            [print(f'A entrada {var} é nula!') for var in self.__dict__.values() if var is None]
            raise ValueError('Há entradas faltando para execução da classe!')

        if len(grid) != dimension:
            raise ValueError('Há discrepância entre a definição da malha e a dimensão adotada para modelagem!')

        if self.dimension != 1:
            raise NotImplementedError('Não está implementada a modelagem de problemas com dimensão diferente de 1! ')

        self.coordinates = self.coordinates.lower()

    def model_parameters(self, eta : float, k : float, phi : float, mu : float, ct : float, lengths : list, area : float, time_list : list, rw = None):
        """
        Define os parâmetros de entrada do modelo.

        eta : float
            É o coeficiente difusivo.
        k : float
            É a permeabilidade.
        phi : float
            É a porosidade.
        mu : float
            É a viscosidade do fluido.
        ct : float
            É a compressibilidade total.
        lengths : [Lx, Ly, Lz] ou [R, h, theta]
            São os comprimentos/tamanhos do modelo. O primeiro comprimento tem preferência sobre os demais.
            Assim, por exemplo, num problema 1D o comprimento Lx será utilizado e Ly e Lz tomados como a área
        well_size : float
            Tamanho do poço. Pode ser um comprimento linear ou radial.
        area : float
            A área pode ser fornecida diretamente, caso não se defina todos os comprimentos necessários.
        time_list : tempos a serem modelados
            Define os tempos em que será feita a modelagem.
        """
        if eta is None: self.eta = k / (phi * mu * ct)
        else: self.eta = eta
        self.k, self.phi, self.mu, self.ct = k, phi, mu, ct
        self.time_list = time_list
        self.rw = rw
        self.p0 = self.ci
        if self.dimension == 1:
            pos = 0
            for cc_type, value in zip(*self.cc):
                if pos == 0:
                    if 'dirichlet' == cc_type.lower():
                        self.pw = value
                    if 'neumann' == cc_type.lower():
                        self.qw = value
                elif pos == 1:
                    if 'dirichlet' == cc_type.lower():
                        self.pe = value
                    if 'neumann' == cc_type.lower():
                        self.qe = value

                else:
                    raise IndexError(
                        'Há discrepância entre os índices dos valores das condições de contorno e seus tipos')
                pos += 1
            if self.coordinates == 'linear':
                self.L = lengths[0]
                if self.rw is None:
                    self.L_list = np.linspace(0.1, lengths[0], self.grid[0])
                else:
                    self.L_list = np.linspace(self.rw, lengths[0], self.grid[0])
                if area is None: self.area = lengths[1] * lengths[2]
                else: self.area = area

            if self.coordinates == 'radial':
                    self.Re = lengths[0]
                    self.h = lengths[1]
                    if self.rw is None:
                        self.R_list = np.linspace(0.1, lengths[0], self.grid[0])
                    else:
                        self.R_list = np.linspace(self.rw, lengths[0], self.grid[0])

    def run(self):
        if self.coordinates == 'linear':
            if self.dimension == 1:
                if all(['dirichlet' == cond.lower() for cond in self.cc[0]]):
                    self.model = 'p_linear_1D_finito'
                    pressures = []
                    for t in self.time_list:

                        pressures.append(p_linear_1D_finito(self.L_list, t, self.pe, self.pw, self.L, self.k, self.phi, self.mu, self.ct))
            if 'neumann' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][1].lower():
                self.model = 'p_linear_1D_realimentacao'
                pressures = []
                for t in self.time_list:
                    pressures.append(
                        p_linear_1D_realimentacao(self.L_list, t, self.p0, self.qw, self.mu, self.L, self.k, self.area, self.phi, self.ct))
            if all(['neumann' == cond.lower() for cond in self.cc[0]]):
                self.model = 'p_linear_1D_selado'
                pressures = []
                for t in self.time_list:
                    pressures.append(
                        p_linear_1D_realimentacao(x=self.L_list, t=t, p0=self.p0, qw=self.qw, mu=self.mu, L=self.L, k=self.k, A=self.area, phi=self.phi, ct=self.ct))
            
        if self.coordinates == 'radial':
            if 'neumann' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][1].lower():
                self.model = 'p_transiente_1D_radial'
                pressures = []
                for t in self.time_list:
                    pressures.append(
                        p_transiente_1D_radial(self.R_list, t, self.p0, self.qw, self.mu, self.h, self.k, self.phi, self.ct))

            if 'dirichlet' == self.cc[0][0].lower() and 'neumann' == self.cc[0][1].lower():
                self.model = 'p_pseudopermanente_1D_radial'
                pressures = []
                for t in self.time_list:
                    pressures.append(
                        p_pseudopermanente_1D_radial(self.R_list, self.Re, self.rw, t, self.p0, self.qe, self.mu, self.h, self.k, self.phi, self.ct))

        self.pressures = pressures

    def postprocess(self, title: str = None, xlim: list = None, times_to_plot: list = None,
                        pos_to_plot: list = None, units: str = 'SI'):
        from PostProcess import plot_p_curves
        from PostProcess import plot_p
        try:
            plot_p_curves(self.L_list, self.pressures, self.time_list, title, xlim, times_to_plot,units)
            plot_p(self.L_list, self.pressures, self.time_list, title, xlim, times_to_plot, pos_to_plot, units)
        except:
            plot_p_curves(self.R_list, self.pressures, self.time_list, title, xlim, times_to_plot,units)
            plot_p(self.R_list, self.pressures, self.time_list, title, xlim, times_to_plot, pos_to_plot, units)


def p_linear_1D_finito(x, t, pe, pw, L, k, phi, mu, ct, N=100):
    eta = k / (phi * mu * ct)
    soma = np.zeros_like(x, dtype=float)
    
    for n in range(1, N + 1):
        a = (n * np.pi) / L
        termo = (np.exp(-(a**2)*eta*t)/n)*np.sin(a*x)
        soma += termo

    p = (pe - pw)*((x/L) + (2/np.pi)*soma) + pw
    
    return p

def p_linear_1D_realimentacao(x, t, p0, qw, mu, L, k, A, phi, ct):
    # INFINITO
    eta = k / (phi * mu * ct)
    a = (qw * mu * L)/(k*A)
    b = 4*eta*t

    p = p0 - a*(np.sqrt(b/(np.pi*L**2))*np.exp(-(x**2)/b) - ((x/L)*erfc(x/np.sqrt(b))))
    
    return p

def p_linear_1D_realimentacao_finito(x, t, p0, qw, mu, L, k, A, phi, ct, N=100):
    # FINITO
    eta = k / (phi * mu * ct)
    soma = np.zeros_like(x, dtype=float)
    for n in range(N):
        lambda_n = (2 * n + 1) * np.pi / (2 * L)
        # Termo da série
        fator_espacial = np.sin(lambda_n * (L - x)) / ((2 * n + 1) ** 2)
        fator_temporal = np.exp(-(lambda_n ** 2) * eta * t)
        soma += fator_espacial * fator_temporal
    termo_estacionario = (qw * mu / (k * A)) * (L - x)
    termo_transitorio = (8 * qw * mu * L) / (k * A * (np.pi ** 2)) * soma
    p = p0 - termo_estacionario + termo_transitorio
    return p

def p_linear_1D_selado(x, t, p0, qw, mu, L, k, A, phi, ct):
    a = (mu*qw)/(k*L*A)
    b = (mu*qw)/(k*A)
    c = qw/(L*A*phi*ct)
    return a*(x**2)/2 - b*x + c*t + p0


def p_transiente_1D_radial(r, t, p0, qw, mu, h, k, phi, ct):
    a = (qw * mu) / (4 * np.pi * k * h)
    b = (phi * mu * ct * r ** 2) / (4 * k * t)
    return p0 + a * expi(-b)


def p_pseudopermanente_1D_radial(r, re, rw, t, p0, qw, mu, h, k, phi, ct):
    a = (qw * mu) / (2 * np.pi * k * h)
    b = 2 * k * t / (phi * mu * ct * re ** 2)
    return p0 - a * (b - np.log(r / rw) + 1 / 2 * (r / re) ** 2 + np.log(re / rw) - 3 / 4)

class Numerical:
    """
    Classe contendo as soluções numéricas.
    Suporta somente problemas 1D

    dimension : int
        Define a quantidade de dimensões espaciais do modelo.

    coordinates : 'Linear' ou 'Radial'
        Configura o tipo de modelagem a ser utilizada.

    ci : Pressão inicial
        Define a condição inicial do problema. É, estritamente, a pressão inicial no domínio do problema.


    cc : [[tipo da CC, tipo da CC],[lado esquerdo, lado direito]]
        São as condições de contorno do problema. Podem ser do tipo Dirichlet (pressão) ou Neumann (vazão).
        Considera-se que o poço está localizado no lado esquerdo e a fronteira do reservatório no lado direito.

    grid : [nx,ny,nz]
        Define o tamanho da malha do problema.

    system_units : 'SI', 'BR', 'USA' ou 'D'
        Informa o sistema de unidades a ser utilizado na modelagem, em que 'D' se refere a sem dimensão.

    """

    def __init__(self, dimension: int, coordinates: str, ci: float, cc: [list, list], grid: list, theta: float,
                 system_units='SI'):
        self.dimension, self.ci, self.cc, self.grid, self.system_units = dimension, ci, cc, grid, system_units
        self.coordinates = coordinates
        self.theta = theta
        if any([None == var for var in self.__dict__.values()]):
            [print(f'A entrada {var} é nula!') for var in self.__dict__.values() if var is None]
            raise ValueError('Há entradas faltando para execução da classe!')

        if len(grid) != dimension:
            raise ValueError('Há discrepância entre a definição da malha e a dimensão adotada para modelagem!')

        if self.dimension != 1 and self.dimension != 2:
            raise NotImplementedError('Não está implementada a modelagem de problemas com dimensão diferente de 1 e 2! ')

        self.coordinates = self.coordinates.lower()

    def model_parameters(self, eta: float, k: float | list, phi: float, mu: float, ct: float, lengths: list, area: float,
                         final_time: float, nt: int, rw=None):
        """
        Define os parâmetros de entrada do modelo.

        eta : float
            É o coeficiente difusivo.
        k : float
            É a permeabilidade.
        phi : float
            É a porosidade.
        mu : float
            É a viscosidade do fluido.
        ct : float
            É a compressibilidade total.
        lengths : [Lx, Ly, Lz] ou [R, h, theta]
            São os comprimentos/tamanhos do modelo. O primeiro comprimento tem preferência sobre os demais.
            Assim, por exemplo, num problema 1D o comprimento Lx será utilizado e Ly e Lz tomados como a área
        well_size : float
            Tamanho do poço. Pode ser um comprimento linear ou radial.
        area : float
            A área pode ser fornecida diretamente, caso não se defina todos os comprimentos necessários.
        final_time : tempo final que a modelagem vai atingir
            Define o tempo em que será feita a modelagem.
        """
        if eta is None and k is float:
            self.eta = k / (phi * mu * ct)
        else:
            self.eta = eta
        self.k, self.phi, self.mu, self.ct = k, phi, mu, ct
        self.final_time = final_time
        self.p0 = self.ci
        self.nt = nt
        self.dt = self.final_time/self.nt
        self.lengths = lengths

        if len(self.grid) == 1:
            self.dx = lengths[0]/self.grid[0]
        if len(self.grid) == 2:
            self.dx = lengths[0] / self.grid[0]
            self.dy = lengths[1] / self.grid[1]

        if self.dimension == 1:
            self.beta = (self.eta * self.dt) / (self.dx ** 2)

            if self.theta == 0:
                if self.beta > 1 / 2:
                    raise ValueError(
                        'Os valores escolhidos para eta, dx e dt não passaram no critério de convergência!')
            if (0 <= self.theta <= 1) is False:
                raise ValueError(f'O valor de theta {self.theta} não pode ser usado!\nEscolha entre 0 e 1!')
            self.rw = rw
            pos = 0
            for cc_type, value in zip(*self.cc):
                if pos == 0:
                    if 'dirichlet' == cc_type.lower():
                        self.pw = value
                    if 'neumann' == cc_type.lower():
                        self.qw = value
                elif pos == 1:
                    if 'dirichlet' == cc_type.lower():
                        self.pe = value
                    if 'neumann' == cc_type.lower():
                        self.qe = value

                else:
                    raise IndexError(
                        'Há discrepância entre os índices dos valores das condições de contorno e seus tipos')
                pos += 1
            if self.coordinates == 'linear':


                if len(self.grid) == 1:
                    self.L = self.lengths[0]
                    self.nx = self.grid[0]
                    if self.rw is None:
                        self.dx = self.L / self.nx
                    else:
                        self.dx = (self.L-self.rw) / self.nx

                # talvez apagar isso
                if self.rw is None:
                    self.L_list = np.linspace(self.dx / 2, self.L - self.dx / 2, self.nx)
                else:
                    self.L_list = np.linspace(self.rw + self.dx / 2, self.L - self.dx / 2, self.nx)
                if area is None:
                    self.area = lengths[1] * lengths[2]
                else:
                    self.area = area

            if self.coordinates == 'radial':
                self.Re = lengths[0]
                self.h = lengths[1]
                if self.rw is None:
                    self.R_list = np.linspace(0.1, lengths[0], self.grid[0])
                else:
                    self.R_list = np.linspace(self.rw, lengths[0], self.grid[0])

            self.Gamma = -self.mu / (self.k * self.area)
        if self.dimension == 2:
            self.beta = 1/(self.phi * self.mu * self.ct)
            if self.theta != 1:
                raise ValueError(f'O valor de theta {self.theta} não pode ser usado!\nEscolha 1!')

            pos = 0
            for cc_type, value in zip(*self.cc):
                if pos == 0:
                    if 'dirichlet' == cc_type.lower():
                        self.pw = value
                    if 'neumann' == cc_type.lower():
                        self.qw = value
                elif pos == 1:
                    if 'dirichlet' == cc_type.lower():
                        self.pe = value
                    if 'neumann' == cc_type.lower():
                        self.qe = value
                elif pos == 2:
                    if 'dirichlet' == cc_type.lower():
                        self.pn = value
                    if 'neumann' == cc_type.lower():
                        self.qn = value
                elif pos == 3:
                    if 'dirichlet' == cc_type.lower():
                        self.ps = value
                    if 'neumann' == cc_type.lower():
                        self.qs = value
                else:
                    raise IndexError(
                        'Há discrepância entre os índices dos valores das condições de contorno e seus tipos')
                pos += 1

            if self.coordinates == 'linear':
                if len(self.grid) == 2:
                    self.Lx = self.lengths[0]
                    self.Ly = self.lengths[1]
                    self.nx = self.grid[0]
                    self.ny = self.grid[1]
                    self.dx = self.Lx / (self.nx - 1)
                    self.dy = self.Ly / (self.ny - 1)
                    self.rx = self.dt / self.dx ** 2
                    self.ry = self.dt / self.dy ** 2
                    self.L_list = [np.linspace(self.dx / 2, self.Lx - self.dx / 2, self.nx),np.linspace(self.dy / 2, self.Ly - self.dy / 2, self.ny)]

        self.time_list = np.linspace(0, self.final_time, self.nt + 1)

    def run(self):
        if self.coordinates == 'linear':
            if self.dimension == 1:
                if self.theta == 0:
                    self.explicita_1D()
                if self.theta != 0:
                    self.implicita_1D()
            if self.dimension == 2:
                if self.theta == 0:
                    pass
                if self.theta != 0:
                    self.implicita_2D()

    def explicita_1D(self):
        self.pressures = np.zeros((self.nt + 1, self.nx))
        for i in range(self.nx):
            self.pressures[0, i] = self.p0

        for n in range(self.nt):
            for i in range(self.nx):
                if i == 0:
                    # Fronteira esquerda
                    if 'dirichlet' == self.cc[0][0].lower():
                        self.pressures[n+1,i] = 4/3 * self.beta *(self.pressures[n,i+1]-3*self.pressures[n,i]+2*self.pw) + self.pressures[n,i]
                    elif 'neumann' == self.cc[0][0].lower():
                        self.pressures[n+1,i] = self.beta * self.pressures[n,i+1] - self.beta * self.pressures[n,i] + self.pressures[n,i] + (self.eta*self.dt*self.Gamma/self.dx) * self.qw
                elif i == self.nx - 1:
                    # Fronteira direita
                    if 'dirichlet' == self.cc[0][1].lower():
                        self.pressures[n+1, i] = 4/3 * self.beta * (self.pressures[n,i-1]-3*self.pressures[n,i]+2*self.pe) + self.pressures[n,i]
                    elif 'neumann' == self.cc[0][1].lower():
                        self.pressures[n+1,i] = self.beta * self.pressures[n,i-1] - self.beta * self.pressures[n,i] + self.pressures[n,i] + (self.eta*self.dt*self.Gamma/self.dx) * self.qe
                else:
                    # Central
                    self.pressures[n + 1, i] = self.beta * self.pressures[n, i + 1] + (1 - 2 * self.beta) * self.pressures[n, i] + self.beta * self.pressures[n, i - 1]

    def implicita_1D(self):
        t = 0
        Nx = 0

        P_matriz_old = np.ones(self.nx)*self.p0 # Pressão inicial do reservatório
        P_matriz = np.zeros((self.nt+1, self.nx))
        P_matriz[Nx,:] = self.p0
        b_matriz = np.zeros(self.nx)

        for n in range(self.nt):
            if n%5 == 0:
                print(f'Progresso: {n/self.nt*100}%')
            Nx = Nx + 1
            A_matriz = np.zeros((self.nx, self.nx))

            for i in range(self.nx):
                # Condições de Contorno:
                if i == 0:  # Fronteira esquerda
                    if 'dirichlet' == self.cc[0][0].lower():
                        A_matriz[i, i] = 1 + 4 * self.beta * self.theta
                        A_matriz[i, i + 1] = - 4 / 3 * self.beta * self.theta
                        b_matriz[i] = (1-4*self.beta*(1-self.theta))*P_matriz_old[i] + 4/3*self.beta*(1-self.theta)*P_matriz_old[i+1]+8/3*self.beta*self.pw

                    elif 'neumann' == self.cc[0][0].lower():
                        A_matriz[i, i] = 1 + self.beta * self.theta
                        A_matriz[i, i + 1] = -self.beta * self.theta
                        b_matriz[i] = (1-self.beta*(1-self.theta))*P_matriz_old[i] + self.beta*(1-self.theta)*P_matriz_old[i+1] + (self.eta*self.dt*self.Gamma/self.dx) * self.qw

                elif i == self.nx - 1:  # Fronteira direita
                    if 'dirichlet' == self.cc[0][1].lower():
                        A_matriz[i, i - 1] = - 4 / 3 * self.beta * self.theta
                        A_matriz[i, i] =1 + 4 * self.beta * self.theta
                        b_matriz[i] = (1-4*self.beta*(1-self.theta))*P_matriz_old[i] + 4/3*self.beta*(1-self.theta)*P_matriz_old[i-1]+8/3*self.beta*self.pe

                    elif 'neumann' == self.cc[0][1].lower():
                        A_matriz[i, i - 1] = -self.beta * self.theta
                        A_matriz[i, i] = 1 + self.beta * self.theta
                        b_matriz[i] = self.beta * (1-self.theta)*P_matriz_old[i-1] + (1-self.beta*(1-self.theta))*P_matriz_old[i] + (self.eta*self.dt*self.Gamma/self.dx) * self.qe

                else:  # Central
                    A_matriz[i, i - 1] = -self.beta*self.theta
                    A_matriz[i, i] = 2*self.beta*self.theta + 1
                    A_matriz[i, i + 1] = -self.beta*self.theta
                    b_matriz[i] = self.beta*(1-self.theta)*P_matriz_old[i+1] +(-2*self.beta*(1-self.theta)+1)*P_matriz_old[i] + self.beta*(1-self.theta)*P_matriz_old[i-1]

            # Qual método vai usar...

            x0 = P_matriz_old # chute inicial
            P_matriz_new = Solve.TDMA(A_matriz, b_matriz)
            P_matriz_old = P_matriz_new.copy()
            P_matriz[n + 1, :] = P_matriz_new
        self.pressures = P_matriz

    def implicita_2D(self):
        t = 0
        P_grid = []
        P_matriz_old = np.ones(self.nx*self.ny)*self.p0 # Pressão inicial do reservatório
        P_matriz = np.zeros((self.nt+1, self.nx*self.ny))
        P_matriz[0,:] = self.p0
        P_grid.append(P_matriz[0].reshape((self.ny, self.nx)))
        b_matriz = np.zeros(self.nx*self.ny)

        for n in range(self.nt):
            if n%5 == 0:
                print(f'Progresso: {n/self.nt*100}%')
            A_matriz = np.zeros((self.nx*self.ny, self.nx*self.ny))

            # Central
            for i in range(1, self.ny - 1):
                for j in range(i * self.nx + 1, (i + 1) * self.nx - 1):
                    k_row = j // self.nx +1
                    k_col = j % self.nx +1
                    k_east = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col + 1])
                    k_west = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col - 1])
                    k_north = 2/(1/self.k[k_row - 1, k_col] + 1/self.k[k_row, k_col])
                    k_south = 2/(1/self.k[k_row + 1, k_col] + 1/self.k[k_row, k_col])
                    A_matriz[j, j] = self.theta*(self.beta*self.rx*(k_east+k_west)+
                                                 self.beta*self.ry*(k_south+k_north)) + 1  # Ap
                    A_matriz[j, j - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                    A_matriz[j, j + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                    A_matriz[j, j - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                    A_matriz[j, j + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                    b_matriz[j] = ((1-self.theta)*self.beta*(self.rx*(k_east*P_matriz_old[j+1]+k_west*P_matriz_old[j-1])
                                    +self.ry*(k_south*P_matriz_old[j+self.nx]+k_north*P_matriz_old[j-self.nx]))+
                                   ((1-self.theta)*(self.beta*(self.rx*(-k_east-k_west)+self.ry*(-k_south-k_north)))+1)*P_matriz_old[j])  # Fonte

            # Fronteira Norte
            for i in range(1, self.nx - 1):
                k_row = 1
                k_col = i+1
                k_east = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col + 1])
                k_west = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col - 1])
                k_north = 2/(1/self.k[k_row - 1, k_col] + 1/self.k[k_row, k_col])
                k_south = 2/(1/self.k[k_row + 1, k_col] + 1/self.k[k_row, k_col])
                if 'dirichlet' == self.cc[0][2].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east+k_west)+self.ry*(3*k_north+k_south)) + 1  # Ap
                    A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                    A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                    A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*(1/3*k_north+k_south) # As
                    b_matriz[i] = (1-self.theta)*(self.beta*(self.rx*(k_east*P_matriz_old[i+1]+k_west*P_matriz_old[i-1])+ self.ry*(1/3*k_north+k_south)*P_matriz_old[self.nx+i]
                                                  +(self.rx*(-k_east-k_west)+self.ry*(-3*k_north-k_south))*P_matriz_old[i]))+P_matriz_old[i]+8/3*self.beta*self.ry*k_north*self.cc[1][2]  # Fonte

                elif 'neumann' == self.cc[0][2].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_west+k_east)+self.ry*(k_south))+1  # Ap
                    A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                    A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                    A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                    b_matriz[i] = (1-self.theta)*(self.beta*(self.rx*(k_east*P_matriz_old[i+1]+k_west*P_matriz_old[i-1])+self.ry*(k_south*P_matriz_old[self.nx+i])
                                                  +(self.rx*(-k_east-k_west)+self.ry*(-k_south)))*P_matriz_old[i])+P_matriz_old[i]  # Fonte

            # Fronteira Sul
            for i in range((self.ny - 1) * self.nx + 1, self.ny * self.nx - 1):
                k_row = i // self.nx + 1
                k_col = i % self.nx + 1
                k_east = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col + 1])
                k_west = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col - 1])
                k_north = 2/(1/self.k[k_row - 1, k_col] + 1/self.k[k_row, k_col])
                k_south = 2/(1/self.k[k_row + 1, k_col] + 1/self.k[k_row, k_col])
                if 'dirichlet' == self.cc[0][3].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east + k_west)+self.ry*(k_north + 3*k_south)) + 1  # Ap
                    A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                    A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                    A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*(k_north + k_south/3)  # An
                    b_matriz[i] = P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*(k_east*P_matriz_old[i+1]+ k_west*P_matriz_old[i-1]-(k_east+k_west)*P_matriz_old[i])+
                                    self.ry*((k_north+k_south/3)*P_matriz_old[i-self.nx]-(k_north+3*k_south)*P_matriz_old[i]))+8/3*self.beta*self.ry*k_south*self.cc[1][3]  # Fonte
                elif 'neumann' == self.cc[0][3].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_west+k_east)+self.ry*k_north) + 1  # Ap
                    A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                    A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                    A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                    b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*(k_east*P_matriz_old[i+1]+k_west*P_matriz_old[i-1]
                                    -(k_east+k_west)*P_matriz_old[i])+self.ry*k_north*(P_matriz_old[i-self.nx]-P_matriz_old[i])))  # Fonte
            # Fronteira Leste
            for i in range(2 * self.nx - 1, (self.ny - 1) * self.nx, self.nx):
                k_row = i // self.nx +1
                k_col = i % self.nx +1
                k_east = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col + 1])
                k_west = 2/(1/self.k[k_row, k_col] + 1/self.k[k_row, k_col - 1])
                k_north = 2/(1/self.k[k_row - 1, k_col] + 1/self.k[k_row, k_col])
                k_south = 2/(1/self.k[k_row + 1, k_col] + 1/self.k[k_row, k_col])
                if 'dirichlet' == self.cc[0][1].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_west+3*k_east)+self.ry*(k_north+k_south)) + 1  # Ap
                    A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*(k_west+k_east/3)  # Aw
                    A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                    A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                    b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_west+k_east/3)*P_matriz_old[i-1]-(k_west+3*k_east)*P_matriz_old[i])
                                    +self.ry*(k_north*P_matriz_old[i-self.nx]+k_south*P_matriz_old[i+self.nx]-(k_north+k_south)*P_matriz_old[i]))+ 8/3*self.beta*self.rx*k_east*self.cc[1][1])  # Fonte
                elif 'neumann' == self.cc[0][1].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*k_west+self.ry*(k_north+k_south)) + 1  # Ap
                    A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                    A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                    A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                    b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_west*(P_matriz_old[i-1]-P_matriz_old[i])+ self.ry*(k_north*P_matriz_old[i-self.nx]
                            +k_south*P_matriz_old[i+self.nx]-(k_north+k_south)*P_matriz_old[i])))  # Fonte

            # Fronteira Oeste
            for i in range(self.nx, (self.ny - 2) * self.nx + 1, self.nx):
                k_row = i // self.nx + 1
                k_col = i % self.nx + 1
                k_east = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col + 1])
                k_west = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col - 1])
                k_north = 2 / (1 / self.k[k_row - 1, k_col] + 1 / self.k[k_row, k_col])
                k_south = 2 / (1 / self.k[k_row + 1, k_col] + 1 / self.k[k_row, k_col])
                if 'dirichlet' == self.cc[0][0].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east+3*k_west)+ self.ry*(k_north+k_south)) + 1  # Ap
                    A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*(k_east+k_west/3)  # Ae
                    A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                    A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                    b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_east+k_west/3)*P_matriz_old[i+1]-(k_east+3*k_west)*P_matriz_old[i])
                        +self.ry*(k_north*P_matriz_old[i-self.nx]+k_south*P_matriz_old[i+self.nx]-(k_north+k_south)*P_matriz_old[i]))+ 8/3*self.beta*self.rx*k_west*self.cc[1][0])  # Fonte
                elif 'neumann' == self.cc[0][0].lower():
                    A_matriz[i, i] = self.theta*self.beta*(self.rx*k_east+self.ry*(k_north+k_south))+1  # Ap
                    A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                    A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                    A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                    b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_east*(P_matriz_old[i+1]-P_matriz_old[i])
                        + self.ry*(k_north*P_matriz_old[i-self.nx]+k_south*P_matriz_old[i+self.nx]-(k_north+k_south)*P_matriz_old[i])))  # Fonte

            # Canto Superior Direito
            i = self.nx - 1
            k_row = i // self.nx + 1
            k_col = i % self.nx + 1
            k_east = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col + 1])
            k_west = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col - 1])
            k_north = 2 / (1 / self.k[k_row - 1, k_col] + 1 / self.k[k_row, k_col])
            k_south = 2 / (1 / self.k[k_row + 1, k_col] + 1 / self.k[k_row, k_col])
            if 'dirichlet' == self.cc[0][1].lower() and 'dirichlet' == self.cc[0][2].lower():
                # Dirichlet em east e dirichlet em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(3*k_east+k_west)+ self.ry*(3*k_north+k_south)) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*(k_west+k_east/3)  # Aw
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*(k_south+k_north/3)  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_west+k_east/3)*P_matriz_old[i-1]-(3*k_east+k_west)*P_matriz_old[i])
                    +self.ry*((k_south+k_north/3)*P_matriz_old[i+self.nx]-(3*k_north+k_south)*P_matriz_old[i]))+8/3*self.beta*(self.rx*k_east*self.cc[1][1]+self.ry*k_north*self.cc[1][2]))  # Fonte
            elif 'neumann' == self.cc[0][1].lower() and 'dirichlet' == self.cc[0][2].lower():
                # Neumann em east e dirichlet em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_west+ self.ry*(3*k_north+k_south)) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*(k_south+k_north/3)  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_west*(P_matriz_old[i-1]-P_matriz_old[i])+self.ry*((k_south+k_north/3)*P_matriz_old[i+self.nx]
                        -(3*k_north+k_south)*P_matriz_old[i]))+8/3*self.beta*self.ry*k_north*self.cc[1][2])  # Fonte
            elif 'dirichlet' == self.cc[0][1].lower() and 'neumann' == self.cc[0][2].lower():
                # Dirichlet em east e neumann em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(3*k_east+k_west)+self.ry*k_south) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*(k_west+k_east/3)  # Aw
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_west+k_east/3)*P_matriz_old[i-1]
                        -(3*k_east+k_west)*P_matriz_old[i])+self.ry*k_south*(P_matriz_old[i+self.nx]-P_matriz_old[i]))+ 8/3*self.beta*self.rx*k_east*self.cc[1][1])  # Fonte
            elif 'neumann' == self.cc[0][1].lower() and 'neumann' == self.cc[0][2].lower():
                # Neumann em east e neumann em north
                A_matriz[i, i] = self.theta*self.beta*(
                    self.rx*k_west + self.ry*k_south
                ) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_west*(P_matriz_old[i-1]-P_matriz_old[i])
                    +self.ry*k_south*(P_matriz_old[i+self.nx]-P_matriz_old[i])))  # Fonte

            # Canto Superior Esquerdo
            i = 0
            k_row = i // self.nx + 1
            k_col = i % self.nx + 1
            k_east = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col + 1])
            k_west = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col - 1])
            k_north = 2 / (1 / self.k[k_row - 1, k_col] + 1 / self.k[k_row, k_col])
            k_south = 2 / (1 / self.k[k_row + 1, k_col] + 1 / self.k[k_row, k_col])
            if 'dirichlet' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][2].lower():
                # Dirichlet em west e dirichlet em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east+3*k_west)+self.ry*(3*k_north+k_south)) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*(k_east+k_west/3)  # Ae
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*(k_south+k_north/3)  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_east+k_west/3)*P_matriz_old[i+1]-(k_east+3*k_west)*P_matriz_old[i])
                    +self.ry*((k_south+k_north/3)*P_matriz_old[i+self.nx]-(3*k_north+k_south)*P_matriz_old[i]))+ 8/3*self.beta*(self.rx*k_west*self.cc[1][0]+ self.ry*k_north*self.cc[1][2]))  # Fonte
            elif 'neumann' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][2].lower():
                # Neumann em west e dirichlet em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_east+self.ry*(3*k_north+k_south)) + 1 # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*(k_south+k_north/3)  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_east*(P_matriz_old[i+1]-P_matriz_old[i])+self.ry*(
                        (k_south+k_north/3)*P_matriz_old[i+self.nx]-(3*k_north+k_south)*P_matriz_old[i]))+8/3*self.beta*self.ry*k_north*self.cc[1][2])  # Fonte
            elif 'dirichlet' == self.cc[0][0].lower() and 'neumann' == self.cc[0][2].lower():
                # Dirichlet em west e neumann em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east+3*k_west)+ self.ry*k_south) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*(k_east+k_west/3)  # Ae
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_east+k_west/3)*P_matriz_old[i+1]-(k_east+3*k_west)*P_matriz_old[i])
                    +self.ry*k_south*(P_matriz_old[i+self.nx]-P_matriz_old[i]))+8/3*self.beta*self.rx*k_west*self.cc[1][0])  # Fonte
            elif 'neumann' == self.cc[0][0].lower() and 'neumann' == self.cc[0][2].lower():
                # Neumann em west e neumann em north
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_east + self.ry*k_south) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                A_matriz[i, i + self.nx] = -self.theta*self.beta*self.ry*k_south  # As
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_east*(P_matriz_old[i+1]-P_matriz_old[i])
                    +self.ry*k_south*(P_matriz_old[i+self.nx]-P_matriz_old[i])))  # Fonte


            # Canto Inferior Direito
            i = self.ny * self.nx - 1
            k_row = i // self.nx + 1
            k_col = i % self.nx + 1
            k_east = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col + 1])
            k_west = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col - 1])
            k_north = 2 / (1 / self.k[k_row - 1, k_col] + 1 / self.k[k_row, k_col])
            k_south = 2 / (1 / self.k[k_row + 1, k_col] + 1 / self.k[k_row, k_col])
            if 'dirichlet' == self.cc[0][1].lower() and 'dirichlet' == self.cc[0][3].lower():
                # Dirichlet em east e dirichlet em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_west+3*k_east)+self.ry*(k_north+3*k_south)) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*(k_west+k_east/3)  # Aw
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*(k_north+k_south/3)  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_west+k_east/3)*P_matriz_old[i-1]-(k_west+3*k_east)*P_matriz_old[i])
                    +self.ry*((k_north+k_south/3)*P_matriz_old[i-self.nx]-(k_north+3*k_south)*P_matriz_old[i]))+ 8/3*self.beta*(self.rx*k_east*self.cc[1][1]+self.ry*k_south*self.cc[1][3]))  # Fonte
            elif 'neumann' == self.cc[0][1].lower() and 'dirichlet' == self.cc[0][3].lower():
                # Neumann em east e dirichlet em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_west+self.ry*(k_north+3*k_south)) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*(k_north+k_south/3)  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_west*(P_matriz_old[i-1]-P_matriz_old[i])
                    + self.ry*((k_north+k_south/3)*P_matriz_old[i-self.nx]-(k_north+3*k_south)*P_matriz_old[i]))+8/3*self.beta*self.ry*k_south*self.cc[1][3])  # Fonte
            elif 'dirichlet' == self.cc[0][1].lower() and 'neumann' == self.cc[0][3].lower():
                # Dirichlet em east e neumann em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_west+3*k_east)+self.ry*k_north) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*(k_west+k_east/3)  # Aw
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_west+k_east/3)*P_matriz_old[i-1]-(k_west+3*k_east)*P_matriz_old[i])
                    + self.ry*k_north*(P_matriz_old[i-self.nx]-P_matriz_old[i]))+8/3*self.beta*self.rx*k_east*self.cc[1][1])  # Fonte
            elif 'neumann' == self.cc[0][1].lower() and 'neumann' == self.cc[0][3].lower():
                # Neumann em east e neumann em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_west+self.ry*k_north) + 1  # Ap
                A_matriz[i, i - 1] = -self.theta*self.beta*self.rx*k_west  # Aw
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_west*(P_matriz_old[i-1]-P_matriz_old[i])
                    + self.ry*k_north*(P_matriz_old[i-self.nx]-P_matriz_old[i])))  # Fonte

            # Canto Inferior Esquerdo
            i = (self.ny - 1) * self.nx
            k_row = i // self.nx + 1
            k_col = i % self.nx + 1
            k_east = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col + 1])
            k_west = 2 / (1 / self.k[k_row, k_col] + 1 / self.k[k_row, k_col - 1])
            k_north = 2 / (1 / self.k[k_row - 1, k_col] + 1 / self.k[k_row, k_col])
            k_south = 2 / (1 / self.k[k_row + 1, k_col] + 1 / self.k[k_row, k_col])
            if 'dirichlet' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][3].lower():
                # Dirichlet em west e dirichlet em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east+3*k_west)+self.ry*(k_north+3*k_south)) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*(k_east+k_west/3)  # Ae
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*(k_north+k_south/3)  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_east+k_west/3)*P_matriz_old[i+1]-(k_east+3*k_west)*P_matriz_old[i])
                    + self.ry*((k_north+k_south/3)*P_matriz_old[i-self.nx]-(k_north+3*k_south)*P_matriz_old[i]))+ 8/3*self.beta*(self.rx*k_west*self.cc[1][0]+ self.ry*k_south*self.cc[1][3]))  # Fonte
            if 'neumann' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][3].lower():
                # Neumann em west e dirichlet em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_east+self.ry*(k_north+3*k_south)) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*(k_north+k_south/3)  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_east*(P_matriz_old[i+1]-P_matriz_old[i])
                    + self.ry*((k_north+k_south/3)*P_matriz_old[i-self.nx]-(k_north+3*k_south)*P_matriz_old[i]))+8/3*self.beta*self.ry*k_south*self.cc[1][3])  # Fonte
            if 'dirichlet' == self.cc[0][0].lower() and 'neumann' == self.cc[0][3].lower():
                # Dirichlet em west e neumann em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*(k_east+3*k_west)+ self.ry*k_north) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*(k_east+k_west/3)  # Ae
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*((k_east+k_west/3)*P_matriz_old[i+1]-(k_east+3*k_west)*P_matriz_old[i])
                    + self.ry*k_north*(P_matriz_old[i-self.nx]-P_matriz_old[i]))+8/3*self.beta*self.rx*k_west*self.cc[1][0])  # Fonte
            if 'neumann' == self.cc[0][0].lower() and 'neumann' == self.cc[0][3].lower():
                # Neumann em west e neumann em south
                A_matriz[i, i] = self.theta*self.beta*(self.rx*k_east + self.ry*k_north) + 1  # Ap
                A_matriz[i, i + 1] = -self.theta*self.beta*self.rx*k_east  # Ae
                A_matriz[i, i - self.nx] = -self.theta*self.beta*self.ry*k_north  # An
                b_matriz[i] = (P_matriz_old[i]+(1-self.theta)*self.beta*(self.rx*k_east*(P_matriz_old[i+1]-P_matriz_old[i])
                    + self.ry*k_north*(P_matriz_old[i-self.nx]-P_matriz_old[i])))  # Fonte

            x0 = P_matriz_old # chute inicial
            P_matriz_new, _, _,_ = Solve(n = 8).Gauss_Seidel(A=A_matriz, B=b_matriz, x0=x0)


            P_matriz_old = P_matriz_new.copy()
            P_matriz[n + 1, :] = P_matriz_new
            P_grid.append(P_matriz[n + 1].reshape((self.ny, self.nx)))
        self.pressures = P_grid

    def compute_error(self, analytical, times_to_compute_err : list = None):
        if times_to_compute_err == None:
            an_pressures = analytical.pressures
            an_times = analytical.time_list
        else:
            all_an_pressures = analytical.pressures
            all_an_times = analytical.time_list
            an_idxs = []
            an_times = []
            for selected_time in times_to_compute_err:
                for idx, time in enumerate(all_an_times):
                    if time == selected_time:
                        an_idxs.append(idx)
                        an_times.append(time)
            an_pressures = np.array([pressure for pressure, time in zip(all_an_pressures, all_an_times) if time in times_to_compute_err])

        if len(self.pressures[0]) == len(an_pressures[0]) is False:
            raise IndexError('Há inconformidade entre os tamanhos da malha de pressão numérica e analítica!')
        # pegar indice certo da numerica e o grid analitico correspondente calcular erros dai
        numerical_idxs = []
        times_found = []
        for an_time in an_times:
            for idx, time in enumerate(self.time_list):
                if time == an_time:
                    numerical_idxs.append(idx)
                    times_found.append(time)
        if numerical_idxs == []:
            raise AttributeError('Não foi possível encontrar tempos correspondentes entre a malha temporal analítica e numérica!')
        new_an_pressures = np.array([pressure for pressure, time in zip(an_pressures, an_times) if time in times_found])
        new_num_pressures = np.array([self.pressures[i] for i in numerical_idxs])
        local_error_grids = new_an_pressures - new_num_pressures
        E_max_abs = []
        for err in local_error_grids:
            E_max_abs.append(np.max(abs(err)))
        soma = 0
        Err_relative = []
        for an, num in zip(new_an_pressures,new_num_pressures):
            for i in range(len(an)):
                soma += abs((an[i] - num[i]) / an[i])
            Err_relative.append(soma)
            soma = 0
        Err_RMSE = []
        for an, num in zip(new_an_pressures,new_num_pressures):
            for i in range(len(an)):
                soma += abs(an[i] - num[i])**2
            Err_RMSE.append((soma/len(an))**(1/2))
            soma = 0
        self.Err_relative = Err_relative
        self.Err_RMSE = Err_RMSE
        self.E_max_abs = E_max_abs
        self.local_error_grids = local_error_grids
        self.times_found = times_found

    def postprocess(self, title : str = None, xlim : list = None, times_to_plot: list = None, pos_to_plot: list = None, units: str = 'SI', cmap = None):
        from PostProcess import plot_p_curves
        from PostProcess import plot_p
        from PostProcess import plot_p2d_map
        from PostProcess import plot_p2d_3Dsurface
        from PostProcess import maps_plot_2D
        if self.dimension == 1:
            try:
                plot_p_curves(self.L_list, self.pressures, self.time_list, title, xlim, times_to_plot,units)
                plot_p(self.L_list, self.pressures, self.time_list, title, xlim, times_to_plot, pos_to_plot, units)
            except:
                plot_p_curves(self.R_list, self.pressures, self.time_list, title, xlim, times_to_plot,units)
                plot_p(self.R_list, self.pressures, self.time_list, title, xlim, times_to_plot, pos_to_plot, units)
        elif self.dimension == 2:
            plot_p2d_map(self.L_list, self.pressures, self.time_list, title, times_to_plot, units)
            plot_p2d_3Dsurface(self.L_list, self.pressures, self.time_list, title, times_to_plot, units)
            maps_plot_2D(self.time_list, times_to_plot, self.L_list, self.pressures, cmap=cmap, units = units)
########################################################################################
def mapa_permeabilidade(nx, ny, k_low, k_up,file_path='.', file_name='mapa_permeabilidade.csv', plot=True):
    blocks = nx * ny
    sample = qmc.LatinHypercube(d=1)
    sample_permiability = sample.random(n=blocks)
    sample_permiability = qmc.scale(sample_permiability, [k_low], [k_up])


    k_grid = sample_permiability.reshape((ny, nx))

    os.makedirs(file_path, exist_ok=True)
    csv_path = os.path.join(file_path, file_name)

    df_permeabilidade = pd.DataFrame(k_grid)
    df_permeabilidade.to_csv(csv_path, index=False, header=False)

    k_fig = None

    if plot:
        k_fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(k_grid, cmap='jet', aspect='auto', origin='lower', extent=[0, nx, 0, ny])

        barra_cores = k_fig.colorbar(im, ax=ax)
        barra_cores.set_label('Permeabilidade (D)')

        ax.set_title('Distribuição de Permeabilidades')
        ax.set_xlabel('Blocos X')
        ax.set_ylabel('Blocos Y')

        f_name, _ = os.path.splitext(file_name)
        fig_path = os.path.join(file_path, f"{f_name}.png")

        k_fig.savefig(fig_path, dpi=300, bbox_inches='tight')

        plt.show()

    return k_grid, k_fig