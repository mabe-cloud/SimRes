##############
#### MAIN ####
##############
import numpy as np
import pandas as pd

import EDH
from EDH import Analytical
from EDH import Numerical
import PostProcess as p

# Dados do problema
p0 = 300 #kgf/cm^2
pe = p0 #kgf/cm^2
pw = p0
pn = p0
ps = p0
B = 1.2 # m^3 / m^3 std
Lx = 1500 #m
Ly = 1500 #m
h = 30 # m
phi = 0.18
mu = 0.8 #cP
ct = 150e-6 #(kgf/cm^2)^-1
NX = 40
NY = 40

# k_grid, _ = EDH.mapa_permeabilidade(NX+2, NY+2, k_low=0.300, k_up=0.5000001) Geracao da malha de permeabilidades

#######################################
######### Validacao do modelo #########
#######################################

k_grid = np.array(pd.read_csv('mapa_permeabilidade.csv', sep=',', header=None))
t_list = [50]
NT = 30
pw = 100
pe = 300

p2_implicita_validacao = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet','Neumann', 'Neumann'],[pw, pe, pn, ps]], grid=[NX, NY],theta=1, system_units='BR')
p2_implicita_validacao.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=None)
# p2_implicita_validacao.run()

t_list_to_plot = np.linspace(0, t_list[-1], NT + 1)
t_list_to_plot = [0, 10, 20, 30, 40, 50]

# p2_implicita_validacao.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR')

####################################################
######### Reservatório com 1 poco produtor #########
####################################################

pocos = [{'tipo': 'produtor', 'pw': 100, 'rw': 0.08,'indice': 740}] # produtor no meio do modelo

k_grid = np.array(pd.read_csv('mapa_permeabilidade.csv', sep=',', header=None))
t_list = [120]
NT = 60
pw = p0
pe = p0

p2_implicita_1pocoprod = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[pw, pe, pn, ps]], grid=[NX, NY],theta=1, system_units='BR')
p2_implicita_1pocoprod.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=pocos)
# p2_implicita_1pocoprod.run()

t_list_to_plot = np.linspace(0, t_list[-1], NT + 1)
t_list_to_plot = [0, 20, 40, 60, 80, 100, 120]

# p2_implicita_1pocoprod.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR')

#############################################################
######### Reservatório com 1 poco produtor e falhas #########
#############################################################

pocos = [{'tipo': 'produtor', 'pw': 100, 'rw': 0.08,'indice': 558}] # produtor na borda do modelo
k_grid = np.array(pd.read_csv('mapa_permeabilidade_falhas.csv', sep=',', header=None))

p2_implicita_falha = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[pw, pe, pn, ps]], grid=[NX, NY],theta=1, system_units='BR')
p2_implicita_falha.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=pocos)
# p2_implicita_falha.run()

t_list_to_plot = np.linspace(0, t_list[-1], NT + 1)
t_list_to_plot = [0, 20, 40, 60, 80, 100, 120]

# p2_implicita_falha.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR')

###################################################################
######### Five spot - 4 pocos injetores e 1 poco produtor #########
###################################################################

pocos = [{'tipo': 'injetor', 'qw': 150, 'indice': 164}, {'tipo': 'injetor', 'qw': 150, 'indice': 196},
         {'tipo': 'injetor', 'qw': 150, 'indice': 1444}, {'tipo': 'injetor', 'qw': 150, 'indice': 1476}, {'tipo': 'produtor', 'pw': 100, 'rw': 0.08,'indice': 740}]

k_grid = np.array(pd.read_csv('mapa_permeabilidade.csv', sep=',', header=None))

p2_implicita_fivespot = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[pw, pe, pn, ps]], grid=[NX, NY],theta=1, system_units='BR')
p2_implicita_fivespot.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=pocos)
p2_implicita_fivespot.run()

t_list_to_plot = np.linspace(0, t_list[-1], NT + 1)
t_list_to_plot = [0, 20, 40, 60, 80, 100, 120]

p2_implicita_fivespot.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR')