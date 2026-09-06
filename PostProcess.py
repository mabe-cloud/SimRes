from matplotlib import pyplot as plt
from pyexpat import model


def post_process(L :list,pressures:list, time_list:list, model : str, title :str):

    plt.figure(figsize=(9, 6))

    for p,t in zip(pressures,time_list):
        plt.plot(L, p, label=f"t = {t} s")

    plt.xlabel('Posição (m)', size=13)
    plt.ylabel("Pressão (Pa)", size=13)
    if title is not None:
        plt.title(title, size=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()