##################
#### NEW MAIN ####
##################
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


tempos = np.linspace(0, 10, 51)

# Caso 1 - Pressão-Pressão

# p_an_1 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Dirichlet', 'Dirichlet'], [pw,pe]], grid=[100], system_units='BR')
# p_an_1.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
# p_an_1.run()
#
# p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw,pe]], grid=[100],theta=0, system_units='BR')
# p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
# p1_explicita.run()
# p1_explicita.compute_error(p_an_1, times_to_compute_err=[10])
# t_list = [1, 5, 10]
# x_list = [p1_explicita.L_list[0], p1_explicita.L_list[int(1/4*100)], p1_explicita.L_list[int(2/4*100)], p1_explicita.L_list[int(3/4*100)], p1_explicita.L_list[99]]
#
# p1_explicita.postprocess(title='Solução Numérica Explícita',times_to_plot=t_list, pos_to_plot=x_list, units='BR')
#
#
# p1_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw,pe]], grid=[100],theta=1, system_units='BR')
# p1_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
# p1_implicita.run()
# p1_implicita.compute_error(p_an_1)
# x_list = [p1_implicita.L_list[0], p1_implicita.L_list[int(1/4*100)], p1_implicita.L_list[int(2/4*100)], p1_implicita.L_list[int(3/4*100)], p1_implicita.L_list[99]]
#
# p1_implicita.postprocess(title='Solução Numérica Implícita',times_to_plot=[1, 5, 10], pos_to_plot=x_list, units='BR')
#
#
# p.maps_plot(time_list=p1_implicita.time_list, t_selected=t_list, x_pos=p1_implicita.L_list,p_an=p_an_1.pressures, p_an_times=p_an_1.time_list, p_explicit=p1_explicita.pressures, p_implicit=p1_implicita.pressures)
#
# p.plot_p(L = p1_implicita.L_list, Pressures = p1_implicita.pressures, time_list=p1_implicita.time_list, title = 'Método Implícito', times_to_plot = t_list, pos_to_plot = x_list, units = 'BR')


# Cálculo de erros
# NXs = [10, 20, 30, 40]
# NT = 10000
# for Nx in NXs:
#     p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
#                              grid=[Nx], theta=0, system_units='SI')
#     p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=10, nt=NT)
#     p1_explicita.run()
#     p1_explicita.compute_error(p_an_1, times_to_compute_err=[10])

# Caso 2 - Pressão/Vazão:

qw = qw_std*B
# p_an_2 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'Dirichlet'], [qw,pe]], grid=[100], system_units='BR')
# p_an_2.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
# p_an_2.run()
# t_list = [1, 5, 10]
# x_list = [p_an_2.L_list[0], p_an_2.L_list[int(1/4*100)], p_an_2.L_list[int(2/4*100)], p_an_2.L_list[int(3/4*100)], p_an_2.L_list[99]]
#
# p_an_2.postprocess(title='Solução Analítica',times_to_plot=t_list, pos_to_plot=x_list, units='BR')
#
#
# p2_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'],[qw,pe]], grid=[100],theta=0, system_units='BR')
# p2_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
# p2_explicita.run()
# p2_explicita.compute_error(p_an_2, times_to_compute_err=[10])
#
# x_list = [p2_explicita.L_list[0], p2_explicita.L_list[int(1/4*100)], p2_explicita.L_list[int(2/4*100)], p2_explicita.L_list[int(3/4*100)], p2_explicita.L_list[99]]
#
# p2_explicita.postprocess(title='Solução Numérica Explícita',times_to_plot=t_list, pos_to_plot=x_list, units='BR')
#
#
# p2_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'],[qw,pe]], grid=[100],theta=1, system_units='BR')
# p2_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=10, nt=5000)
# p2_implicita.run()
# p2_implicita.compute_error(p_an_2)
# x_list = [p2_implicita.L_list[0], p2_implicita.L_list[int(1/4*100)], p2_implicita.L_list[int(2/4*100)], p2_implicita.L_list[int(3/4*100)], p2_implicita.L_list[99]]
#
# p2_implicita.postprocess(title='Solução Numérica Implícita',times_to_plot=[1, 5, 10], pos_to_plot=x_list, units='BR')
#
#
# p.maps_plot(time_list=p2_implicita.time_list, t_selected=t_list, x_pos=p2_implicita.L_list,p_an=p_an_2.pressures, p_an_times=p_an_2.time_list, p_explicit=p2_explicita.pressures, p_implicit=p2_implicita.pressures)
#
# p.plot_p(L = p2_implicita.L_list, Pressures = p2_implicita.pressures, time_list=p2_implicita.time_list, title = 'Método Implícito', times_to_plot = t_list, pos_to_plot = x_list, units = 'BR')


# Caso 3 - Vazão//Vazão:
t_list = [0.1, 1, 5, 10, 20, 30, 40, 50,100]
tempos = np.linspace(0, 100, 501)
NX = 50
NT = 5000
p_an_3 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'], [qw_std, q_e]], grid=[NX], system_units='BR')
p_an_3.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
p_an_3.run()

x_list_an = [p_an_3.L_list[0], p_an_3.L_list[int(1/4*NX)], p_an_3.L_list[int(2/4*NX)], p_an_3.L_list[int(3/4*NX)], p_an_3.L_list[NX-1]]
p_an_3.postprocess(title='Solução Analítica',times_to_plot=t_list, pos_to_plot=x_list_an, units='BR')

p3_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'],[qw_std, q_e]], grid=[NX],theta=0, system_units='BR')
p3_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=t_list[-1], nt=NT)
p3_explicita.run()
p3_explicita.compute_error(p_an_3, times_to_compute_err=t_list)
x_list_num = [p3_explicita.L_list[0], p3_explicita.L_list[int(1/4*NX)], p3_explicita.L_list[int(2/4*NX)], p3_explicita.L_list[int(3/4*NX)], p3_explicita.L_list[NX-1]]
# p3_explicita.postprocess(title='Solução Numérica Explícita',times_to_plot=t_list, pos_to_plot=x_list_num, units='BR')
p.plot_p(L = p3_explicita.L_list, Pressures = p3_explicita.pressures, time_list=p3_explicita.time_list, title = 'Solução Numérica Método Explícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')



p3_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'],[qw_std, q_e]], grid=[NX],theta=1, system_units='BR')
p3_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=t_list[-1], nt=NT)
p3_implicita.run()
p3_implicita.compute_error(p_an_3, times_to_compute_err=t_list)
print('Explícita: ',p3_explicita.Err_RMSE)
print('Implícita: ',p3_implicita.Err_RMSE)
# x_list = [p3_implicita.L_list[0], p3_implicita.L_list[int(1/4*100)], p3_implicita.L_list[int(2/4*100)], p3_implicita.L_list[int(3/4*100)], p3_implicita.L_list[99]]
p.plot_p(L = p3_implicita.L_list, Pressures = p3_implicita.pressures, time_list=p3_implicita.time_list, title = 'Solução Numérica Método Implícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')

p.maps_plot(time_list=p3_implicita.time_list, t_selected=t_list, x_pos=p3_implicita.L_list,p_an=p_an_3.pressures, p_explicit=p3_explicita.pressures, p_implicit=p3_implicita.pressures, p_an_times=p_an_3.time_list)



