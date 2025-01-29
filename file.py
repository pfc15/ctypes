import ctypes
from ctypes import *

M = 4
R = 0.0


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
    for r in range(int(R*10), 10):
        result = _func.sampen2(array, c_int(M), c_double(float(r/10)), c_int(tamanho)) # sampen2(lista, m, r, n)
        list_result = list(cast(result, ctypes.POINTER(ctypes.c_double * M)).contents)
        retorno.append(tuple(list_result))
        print('--'*25)
    return retorno
    
    


def leitura_arquivo(caminho):
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
        numero.append(float(linha[1]))
    print(numero)
    sampen(numero)

def ler_sampentest(caminho):
    lista = []
    with open(caminho, "r") as fp:
        num = fp.readline()
        while len(num)>=1:
            n = float(num)
            lista.append(n)
            num = fp.readline()
            
    return lista

lista = ler_sampentest("./sampentest.txt")
# lista = ler_sampentest("./limpo.txt")


print(sampen(lista))


# SampEn(0,0.2,6782) = 0.000551
# SampEn(1,0.2,6782) = 0.000506
# SampEn(2,0.2,6782) = 0.000501