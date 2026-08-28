import os
import shutil
import numpy as np # type: ignore
from colorama import Fore, Back # type: ignore
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from itertools import combinations_with_replacement
# =========================
def _Info_():
    print("=====")
    print("We have the following methods:")
    methodList = ("Length()", "MKdir()", "Check_out_Word()", "Process_Word()", "ColorList()", "PlotPara()   ")
    for name in methodList:
        print(f"    {name}") 
    print("=====")    
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
        shutil.rmtree(path)
    elif os.path.isdir(path) and not(rm):
        pass
    else:
        os.makedirs(path, exist_ok=True)
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
        pass

    return color

def PlotPara(tickdir='out', labelsize=14, ticksize=12, titlesize=16, nx=False, dx=False):
    plt.rc('xtick', direction=tickdir, labelsize=ticksize, top=True)
    plt.rc('ytick', direction=tickdir, labelsize=ticksize, right=True)
    plt.rc('axes', titlesize=titlesize, labelsize=labelsize)

    if nx:
        plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(nbins=nx))
    else:
        pass
    
    if dx:
        plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(dx))
    else:
        pass

def Sum_Combination(n, Ntot):
    combinations = []
    for combo in combinations_with_replacement(range(1, Ntot+1), n):
        if sum(combo) == Ntot:
            combo = list(combo)
            combo.sort(reverse=True)
            combo = tuple(combo)
            combinations.append(combo)
    return combinations
