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
Lx = 7000 #m
Ly = 7000 #m
h = 30 # m
phi = 0.18
mu = 50 #cP
ct = 150e-6 #(kgf/cm^2)^-1
NX = 40
NY = 40

# Geracao da malha de permeabilidades!! k_grid, _ = EDH.mapa_permeabilidade(NX+2, NY+2, k_low=40, k_up=180)

#######################################
######### Validacao do modelo #########
#######################################
# k_grid, _ = EDH.mapa_permeabilidade(NX+2, NY+2, k_low=110, k_up=110.00000001, file_name='mapa_permeabilidade_MDV.csv')
k_grid = np.array(pd.read_csv('mapa_permeabilidade.csv', sep=',', header=None))
t_list = [240]
NT = 60
pw = 100
pe = 300

p2_implicita_validacao = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet','Neumann', 'Neumann'],[pw, pe, 0, 0]], grid=[NX, NY],theta=1, method='mvf')
p2_implicita_validacao.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=None)
p2_implicita_validacao.run()

t_list_to_plot = [0, 24*2, 24*4, 24*6, 24*8, 24*10] # 0, 2, 4, 6, 8 e 10 dias

p2_implicita_validacao.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR', cmap='plasma')

####################################################
######### Reservatório com 1 poco produtor #########
####################################################

pocos = [{'tipo': 'produtor', 'pw': 100, 'rw': 0.08,'indice': 740}] # produtor no meio do modelo

k_grid = np.array(pd.read_csv('mapa_permeabilidade.csv', sep=',', header=None))
t_list = [240]
NT = 60
pw = p0
pe = p0

p2_implicita_1pocoprod = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[0, 0, 0, 0]], grid=[NX, NY],theta=1, method='mvf')
p2_implicita_1pocoprod.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=pocos)
p2_implicita_1pocoprod.run()

t_list_to_plot = [0, 24*2, 24*4, 24*6, 24*8, 24*10] # 0, 2, 4, 6, 8 e 10 dias

p2_implicita_1pocoprod.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR', cmap='RdYlBu')

#############################################################
######### Reservatório com 1 poco produtor e falhas #########
#############################################################

pocos = [{'tipo': 'produtor', 'pw': 100, 'rw': 0.08,'indice': 558}] # produtor na borda do modelo
k_grid = np.array(pd.read_csv('mapa_permeabilidade_falhas.csv', sep=',', header=None))

p2_implicita_falha = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[0, 0, 0, 0]], grid=[NX, NY],theta=1, method='mvf')
p2_implicita_falha.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=pocos)
p2_implicita_falha.run()

t_list_to_plot = [0, 24*2, 24*4, 24*6, 24*8, 24*10] # 0, 2, 4, 6, 8 e 10 dias

p2_implicita_falha.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR', cmap='PiYG')

###################################################################
######### Five spot - 4 pocos injetores e 1 poco produtor #########
###################################################################

pocos = [{'tipo': 'injetor', 'rw': 0.08,'qw': 550, 'indice': 164}, {'tipo': 'injetor', 'rw': 0.08, 'qw': 550, 'indice': 196},
         {'tipo': 'injetor', 'rw': 0.08, 'qw': 550, 'indice': 1444}, {'tipo': 'injetor', 'rw': 0.08, 'qw': 550, 'indice': 1476}, {'tipo': 'produtor', 'pw': 100, 'rw': 0.08,'indice': 740}]

k_grid = np.array(pd.read_csv('mapa_permeabilidade.csv', sep=',', header=None))

p2_implicita_fivespot = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[0, 0, 0, 0]], grid=[NX, NY],theta=1, method='mvf')
p2_implicita_fivespot.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=pocos)
p2_implicita_fivespot.run()

t_list_to_plot = [0, 24*2, 24*4, 24*6, 24*8, 24*10] # 0, 2, 4, 6, 8 e 10 dias
p2_implicita_fivespot.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR', cmap='Spectral')