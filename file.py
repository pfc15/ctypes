import ctypes
from ctypes import *
import pandas as pd
import plotly
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import os
import numpy as np

M = 5
R = 0.00
CAMINHO_DADOS = "/home/pfc15/Documents/2024.2/pibic/mle/dados"






libc = cdll.LoadLibrary("libc.so.6")
_func = ctypes.CDLL('/home/pfc15/Documents/2024.2/pibic/mle/ctypes/binario.so')

_func.media.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_func.media.restype = POINTER(c_double)
_func.sampen2.argtypes = (ctypes.POINTER(ctypes.c_double), c_int, c_double, c_int)

def media(lista_numeros):
    global _func
    
    tamanho = len(lista_numeros)
    array_type = c_double*tamanho
    array = array_type(*lista_numeros)
    
    
    result = _func.media(array, ctypes.c_int(tamanho))
    list_result = list(cast(result, ctypes.POINTER(ctypes.c_double * tamanho)).contents)
    print(list(list_result))

def sampen(lista_numeros):
    global _func
    tamanho = len(lista_numeros)
    array_type = c_double*tamanho
    array = array_type(*lista_numeros)
    retorno = []
    for r in range(int(R*100), 102, 5):
        if r == 0: r=1
        result = _func.sampen2(array, c_int(M), c_double(float(r/100)), c_int(tamanho)) # sampen2(lista, m, r, n)
        list_result = list(cast(result, ctypes.POINTER(ctypes.c_double * M)).contents)
        retorno.append(tuple(list_result))
        print('--'*25)
    return retorno
    
    


def leitura_arquivo(caminho, coluna):
    """
    função para ler o arquivo BTS G-STUDIO File

    argumentos:
    caminho -- caminho do arquivo a ser lido type: str

    return:
    pandas.DataFrame
    """

    with open(caminho, "r") as fp:
        linha = fp.readline()
        dado = []
        for i in range(16):
            next(fp)
        linha = fp.readline().replace(',','.').split("	")
        while len(linha)!=1:
            linha[-1] = linha[-1][:-1]
            dado.append(linha)
            linha = fp.readline().replace(',','.').split("	")
    numero = list()
    for linha in dado:
        numero.append(float(linha[coluna]))
    return numero

def ler_sampentest(caminho):
    lista = []
    with open(caminho, "r") as fp:
        num = fp.readline()
        while len(num)>=1:
            n = float(num)
            lista.append(n)
            num = fp.readline()
            
    return lista

def cria_grafico(df, nome):

    print(df.head(100))

    # min_y = float(df.min().min())
    max_y = float(df.max().max())+0.005
    max_x = float(df.index.max())
    fig = go.Scatter(x=df.index,y=df,mode="lines", name=nome)
    
    return fig, max_x, max_y

def salvar_grafico(self, fig, caminho_output=os.getcwd(), nome="fig"):
        fig.write_image(os.path.join(caminho_output, f"{nome}.png"))
        fig.write_html(os.path.join(caminho_output, f"{nome}.html"))

for i in range(1, 10):
    lista_df = []
    titulos = ["acc General subplots", "Gyro General subplots", "Rotation General subplots", "Rotation General subplots"]
    subplots_titles = [["X axis", "Y axis", "Z axis"], ["X axis", "Y axis", "Z axis"], ["Roll", "Pitch", "Yaw"]]
    if (i-1)%3==0:
        layout = go.Layout(
            template="simple_white",
            title_text=nome[(i-1)//3],
            width=1440,
            height=1080
        )
        fig1 = make_subplots(rows=3, cols=1, vertical_spacing = 0.2, subplot_titles= subplots_titles[(i-1)%3])
        fig2 = make_subplots(rows=3, cols=1, vertical_spacing = 0.2, subplot_titles= subplots_titles[(i-1)%3])
    
    for nome in os.listdir(os.path.abspath(CAMINHO_DADOS)):
        lista = leitura_arquivo(CAMINHO_DADOS+"/"+nome, 1)
        df = pd.DataFrame(sampen(lista))
        df.index = (df.index*5)/100
        lista_df.append(df)
        print('-=-='*25)


    # Combine DataFrames into a 3D NumPy array
    combined_array = np.array([d.values for d in lista_df])

    # Calculate mean and median along the first axis (across DataFrames)
    mean_df = pd.DataFrame(np.mean(combined_array, axis=0), columns=df.columns, index=df.index)
    median_df = pd.DataFrame(np.median(combined_array, axis=0), columns=df.columns, index=df.index)
    var_df = pd.DataFrame(np.var(combined_array, axis=0), columns=df.columns, index=df.index)

    


    fig_median, max_x, max_y = cria_grafico(median_df, titulos[(i-1)//3]+"Median_Sampen")
    
    fig1.add_trace(median_df, row=(i-1)%3, col=1)
    fig1.update_yaxes(range=[0 , max_y])
    fig1.update_xaxes(range=[0, max_x])
    fig1.update_layout(yaxis_title="Sampen", xaxis_title="R", legend_title_text="M")

    fig_var, max_x, max_y = cria_grafico(var_df, titulos[(i-1)//3]+"_Variancia_Sampen")

    fig2.add_trace(var_df, row=(i-1)%3, col=1)
    fig2.update_yaxes(range=[0 , max_y])
    fig2.update_xaxes(range=[0, max_x])
    fig2.update_layout(yaxis_title="Sampen", xaxis_title="R", legend_title_text="M")

    if i%3==2:
        caminho_output= '/home/pfc15/Documents/2024.2/pibic/mle/ctypes/output'
        salvar_grafico(fig1, caminho_output, nome=titulos[(i-1)//3].replace(" ", "_")+"Median_Sampen")
        salvar_grafico(fig2, caminho_output, nome=titulos[(i-1)//3].replace(" ", "_")+"Variance_Sampen")
