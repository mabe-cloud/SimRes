##############
#### MAIN ####
##############
import numpy as np

import EDH
from EDH import Analytical
from EDH import Numerical
import PostProcess as p

# Dados do problema
p0 = 300 #kgf/cm^2
pe = p0 #kgf/cm^2
pw = 100
pn = p0
ps = p0
B = 1.2 # m^3 / m^3 std
Lx = 100 #m
Ly = 50 #m
phi = 0.18
mu = 0.8 #cP
ct = 150e-6 #(kgf/cm^2)^-1
NX = 6
NY = 4
k_grid, _ = EDH.mapa_permeabilidade(NX+2, NY+2, k_low=0.100, k_up=0.400)

t_list = [1]
NT = 50

p2_implicita = Numerical(dimension=2, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet','Neumann', 'Neumann'],[pw, pe, pn, ps]], grid=[NX, NY],theta=1, system_units='BR')
p2_implicita.model_parameters(eta=None, k=k_grid, phi=phi, mu=mu, ct=ct, lengths=[Lx, Ly],area=None, final_time=t_list[-1], nt=NT)
p2_implicita.run()

t_list_to_plot = np.linspace(0, t_list[-1], NT + 1)
t_list_to_plot = [0,0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,1]

p2_implicita.postprocess(title='Solução Numérica 2D',times_to_plot=t_list_to_plot, units='BR')
# p2_implicita.plot_p(L = p1_explicita.L_list, Pressures = p1_explicita.pressures, time_list=p1_explicita.time_list, title = 'Solução Numérica Método Explícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')