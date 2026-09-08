from matplotlib import pyplot as plt
import numpy as np
def post_process(L :list,pressures:list, time_list:list, title :str = None,xlim : list = None, times_to_plot:list = None):
    # todo - mudar nome para algo como pressure_curves
    plt.figure(figsize=(9, 6))
    if isinstance(pressures, list):
        for p,t in zip(pressures,time_list):
            plt.plot(L, p, label=f"t = {t} s")
    else:
        for i in range(len(pressures)):
            if times_to_plot is not None:
                if any(time_list[i] == time for time in times_to_plot):
                    plt.plot(L, pressures[i], label=f"t = {time_list[i]} s")
            else:
                plt.plot(L, pressures[i], label=f"t = {time_list[i]} s")


    plt.xlabel('Posição (m)', size=13)
    plt.ylabel("Pressão (Pa)", size=13)
    plt.xlim(xlim)
    if title is not None:
        plt.title(title, size=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

def maps_plot(time_list, t_selected, x_pos, p_an, p_explicit, p_implicit, cmap=None, units = 'SI'):
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
        t_unidade = 's'
        p_unidade = 'kgf/cm^2'
        m_unidade = 'm'
    else:
        t_unidade = 's'
        p_unidade = 'psi'
        m_unidade = 'ft'

    n_rows = len(t_selected)
    xgrid = x_pos
    if cmap is None:
        mapa = "YlGnBu"
    else:
        mapa = cmap
    fig, axes = plt.subplots(n_rows, 3, figsize=(14, 3 * n_rows), sharex=True, sharey=True)

    methods = [p_an, p_explicit, p_implicit]
    nomes = ['Analítico', 'Explícito', 'Implícito']

    vmin = min(np.min(p_an), np.min(p_explicit), np.min(p_implicit))
    vmax = max(np.max(p_an), np.max(p_explicit), np.max(p_implicit))

    for i, t_val in enumerate(t_selected):
        for j, method_data in enumerate(methods):
            ax = axes[i, j]
            if j == 0:
                row_data = method_data[i]
                data_2d = np.array([row_data])
                cor = ax.imshow(data_2d, aspect='auto', cmap=mapa,
                                vmin=vmin, vmax=vmax, extent=[xgrid.min(), xgrid.max(), 0, 1])
            else:
                for p in range(len(method_data)):
                    if any(time_list[p] == time for time in t_selected):
                        cor = ax.imshow(np.array([method_data[p]]), aspect='auto', cmap=mapa,
                                        vmin=vmin, vmax=vmax, extent=[xgrid.min(), xgrid.max(), 0, 1])


            if i == 0:
                ax.set_title(nomes[j], fontsize=12)

            if j == 0:
                ax.set_ylabel(f't = {t_val}{t_unidade}', fontsize=11)

            ax.set_yticks([])
            if i == len(t_selected) - 1:
                ax.set_xlabel(f"Posição ({m_unidade})")

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
    fig.colorbar(cor, cax=cbar_ax, label=f"Pressão ({p_unidade})")

    plt.show()
