import numpy as np
from scipy.special import erfc
from scipy.special import expi

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

    grid : (nx,ny,nz)
        Define o tamanho da malha do problema.

    system_units : 'SI', 'BR' ou 'USA'
        Informa o sistema de unidades a ser utilizado na modelagem.

    """

    def __init__(self, dimension : int, coordinates : str, ci : float, cc : [list, list], grid : list, system_units : str):
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

    def model_parameters(self, eta : float, lengths : list, area : float, time_list : list, node_L : list):
        """
        Define os parâmetros de entrada do modelo.

        eta : float
            É o coeficiente difusivo.
        lengths : [Lx, Ly, Lz] ou [Ri, Re, LZ]
            São os comprimentos/tamanhos do modelo. O primeiro comprimento tem preferência sobre os demais.
            Assim, por exemplo, num problema 1D o comprimento Lx será utilizado e Ly e Lz tomados como a área
        area : float
            A área pode ser fornecida diretamente, caso não se defina todos os comprimentos necessários.
        time_list : tempos a serem modelados
            Define os tempos em que será feita a modelagem.
        node_L : [nx, ny, nz]
            Número de blocos espaciais, em cada dimensão
        """
        self.eta = eta
        self.time_list = time_list
        if self.coordinates == 'linear':
            if self.dimension == 1:
                self.L = lengths[0]
                self.L_list = np.linspace(0.1, lengths[0], node_L[0])
                if area is None: self.area = lengths[1] * lengths[2]
                else: self.area = area
                pos = 0
                for cc_type, value in zip(*self.cc):
                    if pos == 0:
                        if 'dirichlet' == cc_type.lower():
                            self.pw = value
                        if 'neumann' == cc_type.lower():
                            self.qw = value
                    elif pos ==1:
                        if 'dirichlet' == cc_type.lower():
                            self.pe = value
                        if 'neumann' == cc_type.lower():
                            self.qe = value

                    else: raise IndexError('Há discrepância entre os índices dos valores das condições de contorno e seus tipos')
                    pos += 1

        self.get_model()

    def get_model(self):
        if self.coordinates == 'linear':
            if self.dimension == 1:
                if all(['dirichlet' == cond.lower() for cond in self.cc[0]]):
                    self.model = 'p1_1D'

    def run(self):
        if self.model == 'p1_1D':
            pressures = []
            for t in self.time_list:
                pressures.append(p1_1D(self.L_list, t, self.pe, self.pw, self.eta, self.L))
        from PostProcess import post_process
        post_process(self.L_list,pressures, self.time_list)



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
    def __init__(self):
        pass