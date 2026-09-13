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



# Caso 1 - Pressão-Pressão

t_list = [0.1, 0.3, 0.5, 1, 3, 5, 7, 10]

tempos = np.linspace(0, t_list[-1], 501)
NX = 80
NT = 5000
p_an_1 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]], grid=[NX], system_units='BR')
p_an_1.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
p_an_1.run()

x_list_an = [p_an_1.L_list[0], p_an_1.L_list[int(1/4*NX)], p_an_1.L_list[int(2/4*NX)], p_an_1.L_list[NX-1]]
p_an_1.postprocess(title='Solução Analítica',times_to_plot=t_list, pos_to_plot=x_list_an, units='BR')

p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw, pe]], grid=[NX],theta=0, system_units='BR')
p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=t_list[-1], nt=NT)
p1_explicita.run()
p1_explicita.compute_error(p_an_1, times_to_compute_err=t_list)
x_list_num = [p1_explicita.L_list[0], p1_explicita.L_list[int(1/4*NX)], p1_explicita.L_list[int(2/4*NX)], p1_explicita.L_list[int(3/4*NX)], p1_explicita.L_list[NX-1]]
# p3_explicita.postprocess(title='Solução Numérica Explícita',times_to_plot=t_list, pos_to_plot=x_list_num, units='BR')
p.plot_p(L = p1_explicita.L_list, Pressures = p1_explicita.pressures, time_list=p1_explicita.time_list, title = 'Solução Numérica Método Explícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')

p1_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'],[pw, pe]], grid=[NX],theta=1, system_units='BR')
p1_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=t_list[-1], nt=NT)
p1_implicita.run()
p1_implicita.compute_error(p_an_1, times_to_compute_err=t_list)

print('Explícita: ',p1_explicita.Err_RMSE)
print('Implícita: ',p1_implicita.Err_RMSE)
p.malha_erros(p1_explicita, p1_implicita, t_list, units='BR')

p.plot_p(L = p1_implicita.L_list, Pressures = p1_implicita.pressures, time_list=p1_implicita.time_list, title = 'Solução Numérica Método Implícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')
t_list = [0.1, 1,5]
p.maps_plot(time_list=p1_implicita.time_list, t_selected=t_list, x_pos=p1_implicita.L_list,p_an=p_an_1.pressures, p_explicit=p1_explicita.pressures, p_implicit=p1_implicita.pressures, p_an_times=p_an_1.time_list, units= 'BR')



# Cálculo de erros
final_time = 10
dt = [0.05,0.025,0.025/2,0.025/4, 0.025/8,0.025/16, 0.025/32] # medindo erro em nx=5
# dt = [1,0.5,0.5/2, 0.5/4, 0.5/8, 0.5/16]
# dt = [10,10/2, 10/4, 10/8, 10/16, 10/32, 10/64, 10/128]
# dt = [t*1000 for t in dt]
NTs = [int(final_time/d) for d in dt]
# NXs = [10,20,30,40,50]
# NXs = [50, 40, 30, 20, 10]
# NXs = [5000, 2500, 1000, 500, 250, 100, 50]
# NXs = [640, 320, 160, 80, 40, 20, 10,5,2]
NXs = [512, 256, 128, 64, 32, 16, 8,4,2]
# Dx perto de 1
dx = [L/x for x in NXs]
Erros = {}
for Nt in NTs:
    print('Em Nt = ', Nt)

    tempos = np.linspace(0, final_time, 501)
    temp_erro = []
    for Nx in NXs:
        p_an_1 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
                            grid=[Nx],
                            system_units='BR')
        p_an_1.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, time_list=tempos)
        p_an_1.run()
        print('Erro com Nx = ', Nx)
        # Explicita
        # p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
        #                          grid=[Nx], theta=0, system_units='SI')
        # p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=final_time, nt=Nt)
        # p1_explicita.run()
        # p1_explicita.compute_error(p_an_1, times_to_compute_err=[10])
        # temp_erro.append(p1_explicita.Err_RMSE)

        # Implicita
        p1_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
                                 grid=[Nx], theta=1, system_units='SI')
        p1_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=final_time, nt=Nt)
        p1_implicita.run()
        p1_implicita.compute_error(p_an_1, times_to_compute_err=[5])
        temp_erro.append(p1_implicita.Err_RMSE[0])
        print(p1_implicita.Err_RMSE)
    Erros[Nt] = temp_erro

p.error_plots(Erros,dx, dt)

# Caso 2 - Pressão/Vazão:

qw = qw_std*B

t_list = [0.1, 0.3, 0.5, 1, 3, 5, 7, 10]
tempos = np.linspace(0, t_list[-1], 501)
NX = 80
NT = 5000
p_an_2 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'], [qw, pe]], grid=[NX], system_units='BR')
p_an_2.model_parameters(eta=eta,k=k,phi=phi,mu=mu,ct=ct,lengths=[L], area=A,time_list=tempos)
p_an_2.run()

# mudando o size
# t_pressure = []
# l_list = []
# for i in range(NX):
#      t_pressure.append(p_an_2.pressures[i][0:NX])
# p_an_2.pressures = t_pressure
# p_an_2.L_list = p_an_2.L_list[0:NX]

x_list_an = [p_an_2.L_list[0], p_an_2.L_list[int(1/4*NX)], p_an_2.L_list[int(2/4*NX)], p_an_2.L_list[NX-1]]
p_an_2.postprocess(title='Solução Analítica',times_to_plot=t_list, pos_to_plot=x_list_an, units='BR')

p2_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'], [qw, pe]], grid=[NX],theta=0, system_units='BR')
p2_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=t_list[-1], nt=NT)
p2_explicita.run()
p2_explicita.compute_error(p_an_2, times_to_compute_err=t_list)


x_list_num = [p2_explicita.L_list[0], p2_explicita.L_list[int(1/4*NX)], p2_explicita.L_list[int(2/4*NX)], p2_explicita.L_list[int(3/4*NX)], p2_explicita.L_list[NX-1]]
p.plot_p(L = p2_explicita.L_list, Pressures = p2_explicita.pressures, time_list=p2_explicita.time_list, title = 'Solução Numérica Método Explícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')

p2_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'], [qw, pe]], grid=[NX],theta=1, system_units='BR')
p2_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L],area=A, final_time=t_list[-1], nt=NT)
p2_implicita.run()
p2_implicita.compute_error(p_an_2, times_to_compute_err=t_list)

print('Explícita: ',p2_explicita.Err_RMSE)
print('Implícita: ',p2_implicita.Err_RMSE)

p.plot_p(L = p2_implicita.L_list, Pressures = p2_implicita.pressures, time_list=p2_implicita.time_list, title = 'Solução Numérica Método Implícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')
p.malha_erros(p2_explicita, p2_implicita, units='BR')
t_list = [0.1, 1,5]
p.maps_plot(time_list=p2_implicita.time_list, t_selected=t_list, x_pos=p2_implicita.L_list,p_an=p_an_2.pressures, p_explicit=p2_explicita.pressures, p_implicit=p2_implicita.pressures, p_an_times=p_an_2.time_list, cmap='magma', units = 'BR')


# Cálculo de erros
final_time = 10
dt = [0.05,0.025,0.025/2,0.025/4, 0.025/8,0.025/16, 0.025/32]
# dt = [t*1000 for t in dt]
NTs = [int(final_time/d) for d in dt]
# NXs = [10,20,30,40,50]
# NXs = [50, 40, 30, 20, 10]
# NXs = [5000, 2500, 1000, 500, 250, 100, 50]
# NXs = [640, 320, 160, 80, 40, 20, 10,5,2]
NXs = [512, 256, 128, 64, 32, 16, 8,4,2]
# Dx perto de 1
dx = [L/x for x in NXs]
Erros = {}
for Nt in NTs:
    print('Em Nt = ', Nt)

    tempos = np.linspace(0, final_time, 501)
    temp_erro = []
    for Nx in NXs:
        p_an_2 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['Neumann', 'dirichlet'], [qw, pe]],
                            grid=[Nx], system_units='BR')
        p_an_2.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, time_list=tempos)
        p_an_2.run()
        print('Erro com Nx = ', Nx)
        # Explicita
        # p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
        #                          grid=[Nx], theta=0, system_units='SI')
        # p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=final_time, nt=Nt)
        # p1_explicita.run()
        # p1_explicita.compute_error(p_an_1, times_to_compute_err=[10])
        # temp_erro.append(p1_explicita.Err_RMSE)

        # Implicita
        p2_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
                                 grid=[Nx], theta=1, system_units='SI')
        p2_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=final_time, nt=Nt)
        p2_implicita.run()
        p2_implicita.compute_error(p_an_2, times_to_compute_err=[5])
        temp_erro.append(p2_implicita.Err_RMSE[0])
        print(p2_implicita.Err_RMSE)
    Erros[Nt] = temp_erro

p.error_plots(Erros,dx, dt)


# # Caso 3 - Vazão//Vazão:
# t_list = [0.1, 1, 5, 10, 20, 30, 40, 50,100]
t_list = [0.1, 0.3, 0.5, 1, 3, 5, 7, 10]
tempos = np.linspace(0, t_list[-1], 501)
NX = 80
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
p.malha_erros(p3_explicita, p3_implicita, units='BR')
p.plot_p(L = p3_implicita.L_list, Pressures = p3_implicita.pressures, time_list=p3_implicita.time_list, title = 'Solução Numérica Método Implícito', times_to_plot = t_list, pos_to_plot = x_list_num, units = 'BR')
t_list = [0.1, 1,5]
p.maps_plot(time_list=p3_implicita.time_list, t_selected=t_list, x_pos=p3_implicita.L_list,p_an=p_an_3.pressures, p_explicit=p3_explicita.pressures, p_implicit=p3_implicita.pressures, p_an_times=p_an_3.time_list, cmap='BuPu', units='BR')


# Cálculo de erros
final_time = 10
# dt = [1,0.5,0.25,0.1,0.05,0.025]
dt = [0.05,0.025,0.025/2,0.025/4, 0.025/8,0.025/16, 0.025/32] # medindo erro em nx=5
# dt = [t*1000 for t in dt]
NTs = [int(final_time/d) for d in dt]
# NXs = [10,20,30,40,50]
# NXs = [50, 40, 30, 20, 10]
# NXs = [5000, 2500, 1000, 500, 250, 100, 50]
NXs = [512, 256, 128, 64, 32, 16, 8,4,2]
# Dx perto de 1
dx = [L/x for x in NXs]
Erros = {}
for Nt in NTs:
    print('Em Nt = ', Nt)

    tempos = np.linspace(0, final_time, 501)

    # p_an_1 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]], grid=[Nx],
    #                     system_units='BR')
    # p_an_1.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, time_list=tempos)
    # p_an_1.run()
    temp_erro = []
    for Nx in NXs:
        p_an_3 = Analytical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'], [qw_std, q_e]],
                            grid=[Nx], system_units='BR')
        p_an_3.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, time_list=tempos)
        p_an_3.run()

        print('Erro com Nx = ', Nx)
        # Explicita
        # p1_explicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['dirichlet', 'dirichlet'], [pw, pe]],
        #                          grid=[Nx], theta=0, system_units='SI')
        # p1_explicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=final_time, nt=Nt)
        # p1_explicita.run()
        # p1_explicita.compute_error(p_an_1, times_to_compute_err=[10])
        # temp_erro.append(p1_explicita.Err_RMSE)

        # Implicita
        p3_implicita = Numerical(dimension=1, coordinates='Linear', ci=p0, cc=[['neumann', 'neumann'], [qw_std, q_e]],
                                 grid=[Nx], theta=1, system_units='SI')
        p3_implicita.model_parameters(eta=eta, k=k, phi=phi, mu=mu, ct=ct, lengths=[L], area=A, final_time=final_time, nt=Nt)
        p3_implicita.run()
        p3_implicita.compute_error(p_an_3, times_to_compute_err=[5])
        temp_erro.append(p3_implicita.Err_RMSE[0])
        print(p3_implicita.Err_RMSE)
    Erros[Nt] = temp_erro

p.error_plots(Erros,dx, dt)

