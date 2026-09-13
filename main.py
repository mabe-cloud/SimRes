import numpy as np
from EDH import Analytical
from EDH import Numerical
import PostProcess as p


# Dados do problema

p0 = 300 #kgf/cm^2
pe = p0 #kgf/cm^2
pw = 150 #kgf/cm^2
qw_std = 400 # m^3 std / dia
q_e = 0
B = 1.2 # m^3 / m^3 std
L = 10000 #m
w = 200 #m
h = 20 #m
A = w*h #m^2
k = 20 #mD 
phi = 0.18
mu = 0.8 #cP
ct = 150e-6 #(kgf/cm^2)^-1
eta = k / (phi * mu * ct)


tempos = np.linspace(0, 10, 50)

# Caso 1 - Pressão-Pressão

p_an_1 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Dirichlet', 'Dirichlet'], [pw,pe]], grid=[100], system_units='SI')
p_an_1.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
p_an_1.run()

p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw,pe]], grid=[100],theta=0, system_units='SI')
p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
p1_explicita.run()


p1_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw,pe]], grid=[100],theta=1, system_units='SI')
p1_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
p1_implicita.run()

t_list = [1, 5, 10]
x_list = [p1_implicita.L_list[0], p1_implicita.L_list[int(1/4*100)], p1_implicita.L_list[int(2/4*100)], p1_implicita.L_list[int(3/4*100)], p1_implicita.L_list[99]]

p.maps_plot(time_list=p1_implicita.time_list, t_selected=t_list, x_pos=p1_implicita.L_list,p_an=p_an_1.pressures, p_explicit=p1_explicita.pressures, p_implicit=p1_implicita.pressures)
p.plot_p(L = p1_implicita.L_list, Pressures = p1_implicita.pressures, time_list=p1_implicita.time_list, title = 'Método Implícito', times_to_plot = t_list, pos_to_plot = x_list, units = 'SI')


# Cálculo de erros

# CAso 2 - Pressão/Vazão:

# Caso 3 - Vazão//Vazão:

p_an_3 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'], [qw_std, q_e]], grid=[100], system_units='SI')
p_an_3.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
p_an_3.run()

p3_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'],[qw_std, q_e]], grid=[100],theta=0, system_units='SI')
p3_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
p3_explicita.run()


p3_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'],[qw_std, q_e]], grid=[100],theta=1, system_units='SI')
p3_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
p3_implicita.run()

t_list = [1, 5, 10]
x_list = [p3_implicita.L_list[0], p3_implicita.L_list[int(1/4*100)], p3_implicita.L_list[int(2/4*100)], p3_implicita.L_list[int(3/4*100)], p3_implicita.L_list[99]]

p.maps_plot(time_list=p3_implicita.time_list, t_selected=t_list, x_pos=p3_implicita.L_list,p_an=p_an_1.pressures, p_explicit=p1_explicita.pressures, p_implicit=p3_implicita.pressures)
p.plot_p(L = p3_implicita.L_list, Pressures = p3_implicita.pressures, time_list=p3_implicita.time_list, title = 'Método Implícito', times_to_plot = t_list, pos_to_plot = x_list, units = 'SI')

