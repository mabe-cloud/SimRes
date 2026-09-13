from matplotlib import pyplot as plt
import numpy as np


def plot_p_curves(L :list,pressures:list, time_list:list, title :str = None,xlim : list = None, times_to_plot:list = None, units = 'SI'):
    if units == "SI":
        t_unidade = 's'
        p_unidade = 'Pa'
        m_unidade = 'm'
    elif units == 'BR':
        t_unidade = 'h'
        p_unidade = 'kgf/cm$^2$'
        m_unidade = 'm'
    else:
        t_unidade = 'h'
        p_unidade = 'psi'
        m_unidade = 'ft'

    plt.figure(figsize=(9, 6),dpi=120)
    # if isinstance(pressures, list):
    #     for p,t in zip(pressures,time_list):
    #         plt.plot(L, p, label=f"t = {t} {t_unidade}")
    # else:
    numerical_idxs = []
    times_found = []
    for time in times_to_plot:
        for idx, num_time in enumerate(time_list):
            if num_time == time:
                numerical_idxs.append(idx)
                times_found.append(time)
    new_num_pressures = np.array([pressures[i] for i in numerical_idxs])
    for pres, time in zip(new_num_pressures, times_found):
        plt.plot(L, pres, label=f"t = {time} {t_unidade}")
        # for i in range(len(pressures)):
        #     if times_to_plot is not None:
        #         if any(time_list[i] == time for time in times_to_plot):
        #             plt.plot(L, pressures[i], label=f"t = {time_list[i]} s")
        #     else:
        #         plt.plot(L, pressures[i], label=f"t = {time_list[i]} s")


    plt.xlabel(f'Posição ({m_unidade})', size=13)
    plt.ylabel(f"Pressão ({p_unidade})", size=13)
    plt.xlim(xlim)
    if title is not None:
        plt.title(title, size=16)
    plt.grid(True, alpha=1)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_p(L: list, Pressures: list, time_list: list, title : str = None, xlim : list = None, times_to_plot: list = None, pos_to_plot: list = None, units = 'SI'):
    """
    Plota gráficos com cortes de p em tempos e posições específicas.
    """
    if units == "SI":
        t_unidade = 's'
        p_unidade = 'Pa'
        m_unidade = 'm'
    elif units == 'BR':
        t_unidade = 'h'
        p_unidade = 'kgf/cm$^2$'
        m_unidade = 'm'
    else:
        t_unidade = 'h'
        p_unidade = 'psi'
        m_unidade = 'ft'

    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=120)
    fig.suptitle(f'{title}',size=16)


    # if times_to_plot is None:
    #     for p,t in zip(Pressures, time_list):
    #         ax1.plot(L, p, label=f"t = {t}{t_unidade}")
    # elif isinstance(Pressures, list):
    #     for p,t in zip(Pressures,time_list):
    #         ax1.plot(L, p, label=f"t = {t}{t_unidade}")
    # else:
    time_numerical_idxs = []
    times_found = []
    for time in times_to_plot:
        for idx, num_time in enumerate(time_list):
            if num_time == time:
                time_numerical_idxs.append(idx)
                times_found.append(time)
    new_num_pressures = np.array([Pressures[i] for i in time_numerical_idxs])
    for pres, time in zip(new_num_pressures, times_found):
        ax1.plot(L, pres, label=f"t = {time} {t_unidade}")
        # for i in range(len(Pressures)):
        #     if any(time_list[i] == time for time in times_to_plot):
        #         ax1.plot(L, Pressures[i], label=f"t = {time_list[i]}{t_unidade}")
    ax1.legend()
    ax1.grid(visible=True, axis='both')
    ax1.set_title(f'Cortes temporais\ndo reservatório')
    ax1.set_xlabel(f'Posição ({m_unidade})', size=14)
    ax1.set_ylabel(f'Pressão ({p_unidade})', size=14)
    if pos_to_plot is None:
        for i,x in enumerate(L):
            p_to_plot = []
            for j in range(len(time_list)):
                p_to_plot.append(Pressures[j][i])
            ax2.plot(time_list, p_to_plot, label=f"pos = {int(x)}{m_unidade}")
    else:
        pos_numerical_idxs = []
        pos_found = []
        for pos in pos_to_plot:
            for idx, num_pos in enumerate(L):
                if num_pos == pos:
                    pos_numerical_idxs.append(idx)
                    pos_found.append(pos)
        new_num_pressures = []
        for idx in pos_numerical_idxs:
            temp_pressure = []
            for i in range(len(time_list)):
                temp_pressure.append(Pressures[i][idx])
            new_num_pressures.append(temp_pressure)

        for press, pos in zip(new_num_pressures, pos_found):
           ax2.plot(time_list, press, label=f"pos = {int(pos)}{m_unidade}")
        # for i,x in enumerate(L):
        #     if any(x == pos for pos in pos_to_plot):
        #         p_to_plot = []
        #         for j in range(len(time_list)):
        #             p_to_plot.append(Pressures[j][i])
        #         ax2.plot(time_list, p_to_plot, label=f"pos = {x}{m_unidade}")
    ax2.legend()
    ax2.set_title(f'Cortes espaciais\ndo reservatório')
    ax2.set_xlabel(f'Tempo ({t_unidade})', size=14)
    ax2.grid(visible=True, axis='both')
    plt.show()

def maps_plot(time_list, t_selected, x_pos, p_an, p_an_times, p_explicit, p_implicit, cmap=None, units = 'SI'):
    """
    Plota mapas comparativos para tempos específicos.
    - Linhas: tempos selecionados (t_selected)
    - Colunas: Métodos (Analítico, Explícito, Implícito)
    - Eixo X: Coordenada espacial x
    """

    if units == "SI":
        t_unidade = 's'
        p_unidade = 'Pa'
        m_unidade = 'm'
    elif units == 'BR':
        t_unidade = 'h'
        p_unidade = 'kgf/cm$^2$'
        m_unidade = 'm'
    else:
        t_unidade = 'h'
        p_unidade = 'psi'
        m_unidade = 'ft'

    n_rows = len(t_selected)
    xgrid = x_pos
    if cmap is None:
        mapa = "YlGnBu"
    else:
        mapa = cmap
    fig, axes = plt.subplots(n_rows, 3, figsize=(14, 3 * n_rows), sharex=True, sharey=True)

    # methods = [p_an, p_explicit, p_implicit]
    nomes = ['Analítico', 'Explícito', 'Implícito']

    vmin = min(np.min(p_an), np.min(p_explicit), np.min(p_implicit))
    vmax = max(np.max(p_an), np.max(p_explicit), np.max(p_implicit))

    an_idxs = []
    an_times_found = []
    for selected_time in t_selected:
        for idx, time in enumerate(p_an_times):
            if time == selected_time:
                an_idxs.append(idx)
                an_times_found.append(time)
    an_pressures = np.array([pressure for pressure, time in zip(p_an, p_an_times) if time in t_selected])
    # Pegar numericas certas
    num_idxs = []
    num_times_found = []
    for selected_time in t_selected:
        for idx, time in enumerate(time_list):
            if time == selected_time:
                num_idxs.append(idx)
                num_times_found.append(time)
    imp_pressures = np.array([pressure for pressure, time in zip(p_implicit, time_list) if time in t_selected])
    exp_pressures = np.array([pressure for pressure, time in zip(p_explicit, time_list) if time in t_selected])
    if an_times_found == num_times_found:
        methods = [an_pressures, exp_pressures, imp_pressures]
        for j, method in enumerate(methods):

            for i, time in enumerate(an_times_found):
                ax = axes[i, j]
                cor = ax.imshow(np.array([method[i]]), aspect='auto', cmap=mapa,
                                vmin=vmin, vmax=vmax, extent=[xgrid.min(), xgrid.max(), 0, 1])
                if i == 0:
                    ax.set_title(nomes[j], fontsize=12)

                if j == 0:
                    ax.set_ylabel(f't = {time}{t_unidade}', fontsize=12)

                ax.set_yticks([])
                if i == len(t_selected) - 1:
                    ax.set_xlabel(f"Posição ({m_unidade})", size=12)
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
        fig.colorbar(cor, cax=cbar_ax, label=f"Pressão ({p_unidade})")
        fig.suptitle('Mapa de pressões para diferentes métodos',size=16)
        plt.show()


    # for i, t_val in enumerate(t_selected):
    #     for j, method_data in enumerate(methods):
    #         ax = axes[i, j]
    #         if j == 0:
    #             row_data = method_data[i]
    #             data_2d = np.array([row_data])
    #             cor = ax.imshow(data_2d, aspect='auto', cmap=mapa,
    #                             vmin=vmin, vmax=vmax, extent=[xgrid.min(), xgrid.max(), 0, 1])
    #         else:
    #             for p in range(len(method_data)):
    #                 if any(time_list[p] == time for time in t_selected):
    #                     cor = ax.imshow(np.array([method_data[p]]), aspect='auto', cmap=mapa,
    #                                     vmin=vmin, vmax=vmax, extent=[xgrid.min(), xgrid.max(), 0, 1])
    #
    #
    #         if i == 0:
    #             ax.set_title(nomes[j], fontsize=12)
    #
    #         if j == 0:
    #             ax.set_ylabel(f't = {t_val}{t_unidade}', fontsize=11)
    #
    #         ax.set_yticks([])
    #         if i == len(t_selected) - 1:
    #             ax.set_xlabel(f"Posição ({m_unidade})")
    #
    # fig.subplots_adjust(right=0.85)
    # cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
    # fig.colorbar(cor, cax=cbar_ax, label=f"Pressão ({p_unidade})")
    #
    # plt.show()

def error_plots(Erros, dx, dt):
    plt.figure()
    for key, time in zip(list(Erros.keys()), dt):
        plt.loglog(dx, Erros[key], marker='o', label=f"$\Delta t$ = {time}")

    plt.xlabel(f'Comprimento Característico - $h(m)$', size=13)
    plt.ylabel(f"Erro na Pressão", size=13)
    plt.grid(True, alpha=0.8,which="both", ls="-")
    plt.legend()
    plt.tight_layout()
    plt.title('Impacto do refinamento da malha espacial')
    plt.show()

    plt.figure()
    for i, h in enumerate(dx):
        erros_h = [Erros[t][i] for t in list(Erros.keys())]
        plt.loglog(dt, erros_h, marker='o', label=f"$h$ = {h} m")
    plt.xlabel('Passo de tempo - $\Delta t$(s)', size=13)
    plt.ylabel('Erro na Pressão', size=13)
    plt.grid(True, alpha=0.8, which="both", ls="-")
    plt.legend(title='Malha Espacial ($h$)')
    plt.tight_layout()
    plt.title('Impacto do refinamento da malha temporal')
    plt.show()

def malha_erros(expli,impli,times = None, units='SI'):

    if units == "SI":
        t_unidade = 's'
        p_unidade = 'Pa'
        m_unidade = 'm'
    elif units == 'BR':
        t_unidade = 'h'
        p_unidade = 'kgf/cm$^2$'
        m_unidade = 'm'
    else:
        t_unidade = 'h'
        p_unidade = 'psi'
        m_unidade = 'ft'

    if times is None:
        times = expli.times_found

    cor = ['indianred', 'navy']
    plt.figure(dpi=100)
    plt.title('Malha de Erros', size=15)
    plt.xlabel(f'Tempo ({t_unidade})', size=13)
    plt.ylabel(f'Erro ({p_unidade})', size=13)
    plt.plot(times,expli.Err_RMSE, ls='-',label='RMSE - Método Explícito', color=cor[0])
    plt.plot(times, impli.Err_RMSE, ls='-', label='RMSE - Método Implícito', color=cor[1])
    plt.plot(times,expli.Err_relative, ls='--',label='Relativo - Método Explícito', color=cor[0])
    plt.plot(times, impli.Err_relative, ls='--', label='Relativo - Método Implícito', color=cor[1])
    plt.legend()
    plt.tight_layout()
    plt.grid(True, alpha=0.8)
    plt.show()