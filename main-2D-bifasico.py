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

#####################################
######### Solucao Analitica #########
#####################################

k_grid_x = np.array(pd.read_csv('mapa_permeabilidade_x.csv', sep=',', header=None))
k_grid_y = np.array(pd.read_csv('mapa_permeabilidade_y.csv', sep=',', header=None))
t_list = [240]
NT = 60
pw = 100
pe = 300

soluca_analitica = Analytical(dimension=2, coordinates='Linear', ci=p0, cc=[['Neumann', 'Neumann','Neumann', 'Neumann'],[0, 0, 0, 0]], grid=[NX, NY],system_units='BR')
soluca_analitica.model_parameters(eta=None, k=[k_grid_x,k_grid_y], phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly, h],area=None, final_time=t_list[-1], nt=NT, pocos=None)
soluca_analitica.run()

t_list_to_plot = [0, 24*2, 24*4, 24*6, 24*8, 24*10] # 0, 2, 4, 6, 8 e 10 dias

p2_implicita_validacao.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR', cmap='plasma')

####################################
######### Solucao Numerica #########
####################################