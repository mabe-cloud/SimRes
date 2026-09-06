from matplotlib import pyplot as plt


def post_process(L :list,pressures:list, time_list:list, title :str = None,xlim : list = None, times_to_plot:list = None):

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
