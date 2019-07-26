from matplotlib.pyplot import plot, fill_between, axis, rcParams
from matplotlib import pyplot as plt
from LabIFSC import linearize, Medida, Unidade
from LabIFSC import M as MM
from pandas import DataFrame
from scipy.optimize import curve_fit
from pathlib import Path
from math import pi
from numpy import arange, array
from sys import argv

### SÓ DÁ PRA USAR HEXADECIMAL PRAS CORES AINDA! ###

um = Medida(1)
pardir = Path(__file__).resolve().parents[1]
OPACIDADE_INCERTEZA = '55'


#======================== CONSTANTES ELEMENTARES ========================#
# Source: Nist CODATA
e = Medida((1.60217662e-19, 0.0000000098e-19), 'C')
h = Medida((6.626070040e-34, 0.000000081e-34), 'J s')
c = Medida(299792458, 'm/s')
me = Medida((9.10938356e-31, 0.00000011e-31), 'kg') # electron mass
mu_0 = Medida((4*pi)*1e-7, 'V s/A m')

redo_flag = '-r' in argv
if redo_flag:
    argv.remove('-r')
    
verbose = '-v' in argv
if verbose:
    argv.remove('-v')

show_flag = '-s' in argv
if show_flag:
    argv.remove('-s')

progress_flag = '-p' in argv
if progress_flag:
    argv.remove('-p')

def df_to_xy(df, labels=None):
# Obtém listas xy a partir de DataFrame. Assim,  DataFrames podem ser passados
# assim como argumento: fit(df, ('x_col', 'y_col'))

    if labels:
        xs = df[labels[0]]
        ys = df[labels[1]]
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
        
    else:
        xs, ys = df[df.columns[:2]].values.T
        plt.xlabel(df.columns[0])
        plt.ylabel(df.columns[1])

    return xs, ys


def plot_reta(a, b=Medida(0), latex=False, cor_incerteza=None,
              label_incerteza='Incerteza', sangria=0, *args, **kwargs):
   
    if not isinstance(a, Medida):
        a = Medida(a)
    if not isinstance(b, Medida):
        b = Medida(b)
       
    if latex:
        label = fr'(${a:latex}$) x + (${b:latex}$)'
    else:
        label = f'({a}) x'
        if b.nominal:
            label += f' + ({b})'

    saved_axis = axis()
    xlim = list(saved_axis[:2])

    xlim[0] -= sangria
    xlim[-1] += sangria

    cor_reta = plot(xlim, [a.nominal*x + b.nominal for x in xlim],
                    label=label, *args, **kwargs)[0]._color

    if cor_incerteza == None:
        cor_incerteza = cor_reta + OPACIDADE_INCERTEZA
        
    fill_between(xlim, [(a.nominal+a.incerteza) * x + b.nominal + b.incerteza for x in xlim],
                 [(a.nominal-a.incerteza) * x + b.nominal-b.incerteza for x in xlim],
                 facecolor=cor_incerteza, label=label_incerteza)

    axis(saved_axis)

    return (a, b)


def hex_to_rgb(hex_color):
    return [int(hex_color[i:i+2], 16) for i in range(1, len(hex_color)) if i%2]

def rgb_to_hex(rgb_color):
    return '#'+''.join([f'{str(hex(i))[2:]:0>2}' for i in rgb_color])

def mod_color(hex_color):
    rgb = hex_to_rgb(hex_color)

    # Achar uma f melhor. Sem graça se ff ou 00.
    f = lambda x: int((x*255)**.5)
    
    return rgb_to_hex(f(x) for x in rgb)


def reta_fit(xs, ys=None, imprimir=False, cor_dados=None, cor_reta=None,
             cor_incerteza=None,
             label_incerteza='Incerteza dos mínimos quadrados',
             usar_incerteza_estatistica=True, latex=False,
             plotar_dados=True, label_dados='Dados experimentais',
             mostrar_origem=False, legenda=True,
             sangria=0, *args, **kwargs):
    """
    Calcula mínimos quadrados, plota com incertezas sombreadas e retorna
    a e b como Medidas.
    """

    if isinstance(xs, DataFrame):
        xs, ys = df_to_xy(xs, ys)

    lin = linearize(xs, ys, imprimir=False)
    a, b, dy, da, db = lin.values()

    if usar_incerteza_estatistica:
        if isinstance(a, Medida):
            #implementar:
            #A = Medida((a.nominal, da.nominal), a.unidade)

            A = Medida((a.nominal, da.nominal))
            A.dimensao = a.dimensao
            A.unidades_originais = a.unidades_originais
            
            B = Medida((b.nominal, db.nominal))
            B.dimensao = b.dimensao
            B.unidades_originais = b.unidades_originais

        else:
            A = Medida((a, da))
            B = Medida((b, db))
            
        a, b = A, B

    if plotar_dados:
        cor_dados = plot(xs, ys, 'o', label=label_dados, color=cor_dados, zorder=10)[0]._color

        if cor_reta == None:
            cor_reta = mod_color(cor_dados)

    if mostrar_origem:
        mostre_origem()
        
    plot_reta(a, b, latex, cor_incerteza, label_incerteza, color=cor_reta, sangria=sangria)

    if legenda:
        plt.legend()

    return (a, b)


def ax_fit(xs, ys=None, plotar_dados=True, cor_dados=None,
           cor_reta=None, cor_incerteza=None, imprimir=False,
           label_incerteza='Incerteza dos mínimos quadrados',
           label_dados='Dados experimentais', latex=False,
           mostrar_origem=False, legenda=True,
           sangria=0, *args, **kwargs):

    if isinstance(xs, DataFrame):
        xs, ys = df_to_xy(xs, ys)

    sao_medidas = isinstance(xs[0], Medida)

    if sao_medidas:
        A_sample = ys[0]/xs[0]
        xs = [x.nominal for x in xs]
        ys = [y.nominal for y in ys]
        
    popt, pcov = curve_fit(lambda x, a: a*x, xs, ys)
    A = Medida((*popt, pcov[0][0] **.5))

    if sao_medidas:
    # Atribuir unidade a A.
        A.dimensao = A_sample.dimensao
        A.unidades_originais = A_sample.unidades_originais

    if plotar_dados:
        cor_dados = plot(xs, ys, 'o', label=label_dados, color=cor_dados, zorder=10)[0]._color

        if cor_reta == None:
            cor_reta = mod_color(cor_dados)

    if mostrar_origem:
        mostre_origem()
        
    plot_reta(A, latex=latex, cor_incerteza=cor_incerteza,
              label_incerteza=label_incerteza, color=cor_reta, sangria=sangria)

    if legenda:
        plt.legend()
        
    return A


def alg_sig(s):
    """ Retorna a casa decimal do primeiro algarismo significativo. """
    # kinda cumbersome
    # breaks if s = 0
    # Problemático se o número termina em 0, tipo 3000+/-5

    flag = True
    p = None
    first = None
    
    for i, c in enumerate(s[::-1]):
        if flag and c not in ['0','.']:
            flag = False
            first = i
        elif c == '.':
            p = i

        if p != None and first != None:
            break
    else:
        p = 0

    if p:
        return 10**-p
        
    return 10**(first-p)


def compare(medida1, medida2):
    if medida1 == medida2: return 1
    elif medida1 != medida2: return -1
    else: return 0


def err_stat(ref, query):
    print (f'{"Medida":>20}:', query)
    print (f'{"Referência":>20}:', ref)

    err = abs(query-ref)
    sum_stdev = ref.incerteza + query.incerteza
    
    print (f'{"Erro absoluto":>20}: {err.nominal:.4f}')
    print (f'{"Erro percentual":>20}: {100*err.nominal/query.nominal:.4f}%')
    print(f'{"Erro / soma desvios":>20}: {err.nominal/sum_stdev:.4f}')
    print (f'{"Equivalência":>20}:', ['inconclusiva.', 'confirmada.', 'negada.'][compare(ref, query)])

    
def mostre_origem():
    ax = axis()
    axis([0, ax[1], 0, ax[3]])


def plot_func(func, rang=None, densidade=100):
    ax = axis()
    if rang == None:
        rang = arange(ax[0], ax[1], (ax[1]-ax[0])/densidade)

    plot(rang, func(rang))
    axis(ax)

def M(lista, incerteza=0, unidade=None):
    return array(MM(list(lista), unidade=unidade, incerteza=incerteza))

def converta(lista, unidade):
    return array([i.converta(unidade) for i in lista])

def projetor_look():
    rcParams['savefig.dpi'] = 300
    rcParams['figure.figsize'] = (16, 9)
    rcParams['font.size'] = 22
    rcParams['lines.linewidth'] = 3
    rcParams['lines.markersize'] = 12
    rcParams['axes.grid'] = True
    # print(rcParams.keys())
    return rcParams
