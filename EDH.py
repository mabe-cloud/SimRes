import numpy as np
from scipy.special import erfc
from scipy.special import expi
from Solvers import Solve

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
                    self.model = 'p1_1D'
                    pressures = []
                    for t in self.time_list:
                        pressures.append(p1_1D(self.L_list, t, self.pe, self.pw, self.eta, self.L))
            if 'neumann' == self.cc[0][0].lower() and 'dirichlet' == self.cc[0][1].lower():
                self.model = 'p2_1D'
                pressures = []
                for t in self.time_list:
                    pressures.append(
                        p2_1D(self.L_list, t, self.p0, self.qw, self.mu, self.L, self.k, self.area, self.phi, self.ct))
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
    def postprocess(self, title : str = None, xlim : list = None, times_to_plot: list = None):
        from PostProcess import post_process
        try:
            post_process(self.L_list, self.pressures, self.time_list, title, xlim, times_to_plot)
        except:
            post_process(self.R_list, self.pressures, self.time_list, title, xlim, times_to_plot)



def p1_1D(x, t, pe, pw, eta,L, N=100):
    soma = np.zeros_like(x, dtype=float)

    for n in range(1, N + 1):
        a = (n * np.pi) / L
        termo = (np.exp(-(a ** 2) * eta * t) / n) * np.sin(a * x)
        soma += termo

    p = (pe - pw) * ((x / L) + (2 / np.pi) * soma) + pw

    return p


def p2_1D(x, t, p0, qw, mu, L, k, A, phi, ct):
    eta = k / (phi * mu * ct)
    a = (qw * mu * L) / (k * A)
    b = 4 * eta * t

    p = p0 - a * (np.sqrt(b / (np.pi * L ** 2)) * np.exp(-(x ** 2) / b) - ((x / L) * erfc(x / np.sqrt(b))))

    return p


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

        if self.dimension != 1:
            raise NotImplementedError('Não está implementada a modelagem de problemas com dimensão diferente de 1! ')

        self.coordinates = self.coordinates.lower()

    def model_parameters(self, eta: float, k: float, phi: float, mu: float, ct: float, lengths: list, area: float,
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
        if eta is None:
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
        self.beta = (self.eta * self.dt) / (self.dx ** 2)


        if self.theta == 0:
            if self.beta > 1 / 2:
                raise ValueError('Os valores escolhidos para eta, dx e dt não passaram no critério de convergência!')
        if (0<=self.theta <=1) is False:
            raise ValueError(f'O valor de theta {self.theta} não pode ser usado!\nEscolha entre 0 e 1!')
        self.rw = rw
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

                if len(self.grid) == 1:
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
        self.time_list = np.linspace(0, self.final_time, self.nt + 1)

    def run(self):
        if self.coordinates == 'linear':
            if self.dimension == 1:
                if self.theta == 0:
                    self.explicita_1D()
                if self.theta != 0:
                    self.implicita_1D()
        self.compute_error()
    def matriz(self):
        """
        Faz a montagem da matriz
        """
        pass
        # if self.dimension == 1:
        #     if self.coordinates == 'linear':

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

        # while t < self.final_time:
        #     Nx = Nx + 1
        #     A_matriz = np.zeros((self.nx, self.nx))
        #
        #     for i in range(self.nx):
        #         # Condições de Contorno:
        #         if i == 0:  # Fronteira esquerda
        #             if 'dirichlet' == self.cc[0][0].lower():
        #                 A_matriz[i, i] = 1 + 4 * self.beta * self.theta
        #                 A_matriz[i, i + 1] = - 4 / 3 * self.beta * self.theta
        #                 b_matriz[i] = (1-4*self.beta*(1-self.theta))*P_matriz_old[i] + 4/3*self.beta*(1-self.theta)*P_matriz_old[i+1]+8/3*self.beta*self.pw
        #
        #             elif 'neumann' == self.cc[0][0].lower():
        #                 A_matriz[i, i] = 1 + self.beta * self.theta
        #                 A_matriz[i, i + 1] = -self.beta * self.theta
        #                 b_matriz[i] = (1-self.beta*(1-self.theta))*P_matriz_old[i] + self.beta*(1-self.theta)*P_matriz_old[i+1] + (self.eta*self.dt*self.Gamma/self.dx) * self.qw
        #
        #         elif i == self.nx - 1:  # Fronteira direita
        #             if 'dirichlet' == self.cc[0][1].lower():
        #                 A_matriz[i, i - 1] = - 4 / 3 * self.beta * self.theta
        #                 A_matriz[i, i] =1 + 4 * self.beta * self.theta
        #                 b_matriz[i] = (1-4*self.beta*(1-self.theta))*P_matriz_old[i] + 4/3*self.beta*(1-self.theta)*P_matriz_old[i-1]+8/3*self.beta*self.pe
        #
        #             elif 'neumann' == self.cc[0][1].lower():
        #                 A_matriz[i, i - 1] = -self.beta * self.theta
        #                 A_matriz[i, i] = 1 + self.beta * self.theta
        #                 b_matriz[i] = self.beta * (1-self.theta)*P_matriz_old[i-1] + (1-self.beta*(1-self.theta))*P_matriz_old[i] - (self.eta*self.dt*self.Gamma/self.dx) * self.qe
        #
        #         else:  # Central
        #             A_matriz[i, i - 1] = -self.beta*self.theta
        #             A_matriz[i, i] = 2*self.beta*self.theta + 1
        #             A_matriz[i, i + 1] = -self.beta*self.theta
        #             b_matriz[i] = self.beta*(1-self.theta)*P_matriz_old[i+1] +(-2*self.beta*(1-self.theta)+1)*P_matriz_old[i] + self.beta*(1-self.theta)*P_matriz_old[i-1]
        #
        #     # Qual método vai usar...
        #
        #     x0 = P_matriz_old # chute inicial
        #     P_matriz_new = Solve.TDMA(A_matriz, b_matriz)
        #     P_matriz_old = P_matriz_new.copy()
        #     t = t + self.dt
        #     P_matriz[Nx, :] = P_matriz_new
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

    def compute_error(self):
        try:
            time_list = np.linspace(0, self.final_time, self.nt + 1)
            analitical = Analytical(self.dimension, self.coordinates, self.ci, self.cc, self.grid, self.system_units)
            analitical.model_parameters(self.eta, self.k, self.phi, self.mu, self.ct, self.lengths, self.area,
                         time_list, self.rw)
            analitical.run()
            an_pressures = analitical.pressures
            E_max = 0
            Err_relativo = []
            soma = 0
            for n in range(len(self.pressures)):
                for i in range(len(self.pressures[0])):
                    soma += abs((an_pressures[n][i] - self.pressures[n][i]) / an_pressures[n][i])
                    E_max = max(E_max,abs(an_pressures[n][i] - self.pressures[n][i]))
                Err_relativo.append(soma)
                soma = 0
            Err_RMSE = []
            for n in range(len(self.pressures)):
                for i in range(len(self.pressures[0])):
                    soma += abs(an_pressures[n][i] - self.pressures[n][i])**2
                Err_RMSE.append((soma/len(self.pressures[0]))**(1/2))
                soma = 0




        except:
            print('Não foi encontrada uma solução analítica para computar os erros!')

    def postprocess(self, title : str = None, xlim : list = None, times_to_plot: list = None):
        from PostProcess import post_process
        try:
            post_process(self.L_list, self.pressures, self.time_list, title, xlim, times_to_plot)
        except:
            post_process(self.R_list, self.pressures, self.time_list, title, xlim, times_to_plot)
