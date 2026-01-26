import os
import numpy as np # type: ignore
import re
from colorama import Fore, Back # type: ignore
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement
# =========================
def Length(a:np.array):
    return (np.sqrt(np.sum(a**2)))

def M_Rotate(deg):
    rad = deg/180*np.pi
    return np.array([[np.cos(rad), -np.sin(rad), 0 ],
                     [np.sin(rad),  np.cos(rad), 0 ],
                     [     0     ,       0     , 1 ]])
# =========================
def MKdir(path, rm):
    if os.path.isdir(path) and rm:
        os.system(f"rm -r {path}")
        os.mkdir(path)
    elif os.path.isdir(path) and not(rm):
        pass
    else:
        os.mkdir(path)
# =========================
def Check_out_Word(word: str):
    print(Fore.BLACK, Back.RED + word + Back.RESET, Fore.RESET)

def Process_Word(word: str):
    print(Fore.YELLOW, word , Fore.RESET)
# =========================  
def ColorList(p=False):
    color = ['#00107f', '#7f7f00', '#ff7a00', '#a000c8', '#ff0000', '#545659', '#79317b', 'green', '#373737']
    if p:
        plt.figure(figsize=(10, 1))
        for i, c in enumerate(color):
            plt.scatter(i, 0, color=c)
        plt.grid('--')
        plt.yticks([])
        plt.show()
    else:
        print("Loading color ...")

    return color

def Sum_Combination(n, Ntot):
    combinations = []
    for combo in combinations_with_replacement(range(1, Ntot+1), n):
        if sum(combo) == Ntot:
            combo = list(combo)
            combo.sort(reverse=True)
            combo = tuple(combo)
            combinations.append(combo)
    return combinations
