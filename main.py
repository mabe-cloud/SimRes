##################
#### NEW MAIN ####
##################
import numpy as np
from EDH import Analytical
from EDH import Numerical
import PostProcess as p

p0 = 50000
pe = p0
pw = 1000
L = 1000
e = 0.2
A = L*e
k = 0.01
phi = 0.20
mu = 0.4
ct = 1e-5
eta = k / (phi * mu * ct)


#%% Regime Permanente ---------------------------------------------------------
tempos = [50, 100, 200]
tempos = [1, 2, 3, 4, 5]

Analitica = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Dirichlet', 'Dirichlet'], [pw,pe]], grid=[100], system_units='SI')
Analitica.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
Analitica.run()
Analitica.postprocess(title="Solução Regime Permanente Linear")

tempos = [1, 2, 3, 4, 5]
qw = 80
# Analitica = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'Dirichlet'], [qw,pe]], grid=[500], system_units='SI')
# Analitica.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
# Analitica.run()
# Analitica.postprocess(title="Regime Transiente Linear")

tempos_longos = [10, 15, 20, 25, 30, 40, 50, 70, 100]

# Analitica.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos_longos)
# Analitica.run()
# Analitica.postprocess(title="Regime Pseudopermanente Linear")

Numerico1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw,pe]], grid=[100],theta=0, system_units='SI')
Numerico1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=5, nt=5000)
Numerico1_explicita.run()
Numerico1_explicita.postprocess(title="Explícito - Regime Permanente Linear",times_to_plot=tempos)

Numerico1_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw,pe]], grid=[100],theta=1, system_units='SI')
Numerico1_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=5, nt=5000)
Numerico1_implicita.run()
Numerico1_implicita.postprocess(title="Implícito - Regime Permanente Linear",times_to_plot=tempos)

p.maps_plot(time_list=Numerico1_implicita.time_list, t_selected=tempos, x_pos=Numerico1_implicita.L_list,p_an=Analitica.pressures, p_explicit=Numerico1_explicita.pressures, p_implicit=Numerico1_implicita.pressures)

# Numerico2_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'],[qw,pe]], grid=[2000],theta=0, system_units='SI')
# Numerico2_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L*10],area=A, final_time=100, nt=200000)
# Numerico2_explicita.run()
# Numerico2_explicita.postprocess(title="Explícito - Regime Transiente Linear", times_to_plot=tempos_longos, xlim=[0,1000])

# Numerico2_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'],[qw,pe]], grid=[2000],theta=1, system_units='SI')
# Numerico2_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L*10],area=A, final_time=100, nt=2000)
# Numerico2_implicita.run()
# Numerico2_implicita.postprocess(title="Implícita - Regime Transiente Linear", times_to_plot=tempos_longos, xlim=[0,1000])
#
# # Dados
# pe = 300 * 98066.5 # kgf/cm2 para Pa
# pw = 150 * 98066.5 # kgf/cm2 para Pa
# rw = 1 # m
# re = 500 # m
# k = 20 * 9.869e-16 # md para m2
# phi = 0.18
# mu = 0.8 *1e-3 # cp para Pa.s
# ct = 150e-6 * 1.0197e-5 # (kgf/cm2)^-1 para Pa^-1
# N = 100
# h = 20 # altura da formacao - m
# # Dados para o transiente
# Bo = 1.2 # m3 / m3std
# qw_std = 400 * 1/86400 # m3std/dia para m3std/s
# p0 = pe # pressão inicial
# qw = Bo * qw_std
# tempos = [60*30, 60*60*3, 60*60*12, 86400*2, 4*86400, 8*86400]
# Analitica_radial = Analytical(dimension=1, coordinates='Radial', ci=p0, cc=[['Neumann', 'Dirichlet'], [qw,pe]], grid=[500], system_units='SI')
# Analitica_radial.model_parameters(eta=None, k=k, phi=phi, mu=mu, ct=ct, lengths=[re,h], area=A,time_list=tempos, rw=rw)
# Analitica_radial.run()
# Analitica_radial.postprocess(title='Solução Transiente Radial 1D')
#
# tempos = [60*30, 86400, 86400*3, 7*86400, 14*86400, 30*86400]
# Analitica_radial = Analytical(dimension=1, coordinates='Radial', ci=p0, cc=[['Dirichlet', 'Neumann'], [p0,qw]], grid=[500], system_units='SI')
# Analitica_radial.model_parameters(eta=None, k=k, phi=phi, mu=mu, ct=ct, lengths=[re,h], area=A,time_list=tempos, rw=rw)
# Analitica_radial.run()
# Analitica_radial.postprocess(title='Solução Pseudopermanente Radial 1D')