# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import glob
import numpy as np
import os
import visualization
import xyPlot
import displayGroupOdbToolset as dgo

wd = os.getcwd()


def extractCurves(databasePath):
    odbPath = databasePath
    o1 = session.openOdb(name=odbPath)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    odb = session.odbs[odbPath]
    session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
        NODAL, ((COMPONENT, 'RF3'), )), ('U', NODAL, ((COMPONENT, 'U3'), )), ), 
        nodeSets=("RP-LOAD", ))
    xy1 = session.xyDataObjects['U:U3 PI: ASSEMBLY N: 1']
    xy2 = session.xyDataObjects['RF:RF3 PI: ASSEMBLY N: 1']
    xy3 = combine(xy1, 0.001*xy2)
    xy3.setValues(sourceDescription='combine("U:U3 PI: ASSEMBLY N: 1",0.001*"RF:RF3 PI: ASSEMBLY N: 1")')
    tmpName = xy3.name
    outName = odbPath[len(wd)+3:-4]
    session.xyDataObjects.changeKey(tmpName, outName)
    del session.xyDataObjects['U:U3 PI: ASSEMBLY N: 1']
    del session.xyDataObjects['RF:RF3 PI: ASSEMBLY N: 1']
    print("{0:} exported!".format(outName))






# DATA WRITTING FUNCTIONS ***************************************************


def save_xy_data_to_txt(
    xy_data, 
    output_filename='xy_data_output.txt', 
    directory=None,
    header=None,
    delimiter='\t'):
    try:
        # Define o diretório de saída
        if directory is None:
            directory = os.getcwd()
        
        # Cria o diretório se não existir
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        # Caminho completo do arquivo
        full_path = os.path.join(directory, output_filename)
        
        # Abre o arquivo para escrita
        with open(full_path, 'w') as file:
            # Escreve cabeçalho, se fornecido
            if header:
                file.write('{0}\n'.format(header))
            
            # Escreve dados
            for row in xy_data:
                # Converte cada elemento para string
                row_str = delimiter.join(str(val) for val in row)
                file.write('{0}\n'.format(row_str))
        
        print('Arquivo salvo com sucesso: {0}'.format(full_path))
        return full_path
    
    except IOError as error:
        print('Erro ao salvar arquivo: {0}'.format(error))
        return None
    except Exception as error:
        print('Erro inesperado: {0}'.format(error))
        return None



def exportLoadCurves():
    so = session.xyDataObjects
    sk = so.keys()
    output_Path = os.getcwd()+"/LoadResults"
    if not os.path.exists(output_Path):
            os.makedirs(output_Path)
    for k in sk:
        ls = [[v[0], v[1]] for v in so[k].data]
        save_xy_data_to_txt(ls, 
                            output_filename='LD-{0:}.txt'.format(k.replace("-ts-AN","")), 
                            directory=output_Path,
                            header="Step\tLoad",
                            delimiter='\t')





def read_xy_data_from_txt(
    filepath, 
    delimiter='\t', 
    skip_header=False):
    """
    Lê dados XY de um arquivo de texto.
    
    Parâmetros:
    -----------
    filepath : str
        Caminho completo do arquivo
    delimiter : str, opcional
        Delimitador entre colunas
    skip_header : bool, opcional
        Pula a primeira linha se for cabeçalho
    
    Retorna:
    --------
    List[List[float]]
        Dados XY lidos do arquivo
    """
    try:
        xy_data = []
        
        with open(filepath, 'r') as file:
            # Pula cabeçalho se necessário
            if skip_header:
                next(file)
            
            # Lê e converte linhas
            for line in file:
                # Remove whitespaces e quebras de linha
                line = line.strip()
                
                # Pula linhas vazias
                if not line:
                    continue
                
                # Converte linha para lista de floats
                row = [float(val) for val in line.split(delimiter)]
                xy_data.append(row)
        
        return xy_data
    
    except FileNotFoundError:
        print('Arquivo não encontrado: {0}'.format(filepath))
        return []
    except ValueError as error:
        print('Erro na conversão de dados: {0}'.format(error))
        return []



# Configurações adicionais
def configure_output(base_dir=None):
    if base_dir is None:
        base_dir = os.path.join(os.getcwd(), 'output')
    
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    
    print('Diretório configurado: {0}'.format(base_dir))
    return base_dir

# END OF DATA WRITTING FUNCTIONS ********************************************



# BEGIN OF HASHIN DATA EXTRACTION =========================================================================

def getH(odbPath):
    o1 = session.openOdb(name=odbPath)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF, ))
    session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
    odb = session.odbs[odbPath]
    output_name = "HS-"+odbPath[len(os.getcwd())+2:-4]
    sl  =     session.xyDataListFromField(odb=odb, outputPosition=INTEGRATION_POINT, 
        variable=(('HSNFCCRT', INTEGRATION_POINT), ), elementSets=(
        "REPAIR-1.SET-REPAIR", ))
    output_Path = os.getcwd()+"/HashinOutputs"
    if not os.path.exists(output_Path):
        os.mkdir(output_Path)
    #
    ls = []
    for i in range(len(sl)):
        tmp = [v[1] for v in sl[i].data]
        ls.append(tmp) 
    #
    #
    arr = np.asarray(ls)
    hs = [[j,np.max(arr[:,j])] for j in range(len(arr[0,:]))]
    #
    for keyLabel in session.xyDataObjects.keys():
        if "HSN" in keyLabel:
            del session.xyDataObjects[keyLabel]
    #
    save_xy_data_to_txt(
    hs, 
    output_filename='{0:}.txt'.format(output_name), 
    directory=output_Path,
    header="Step\tTime",
    delimiter='\t')



def getMax():
    so = session.xyDataObjects
    sk = so.keys()
    lst = [[sk[i], max([v[1] for v in so[sk[i]].data])] for i in range(len(sk))]
    for l in lst: 
        print("{0:}\t{1:0.4f}".format(l[0], l[1]))



def getHashin():
    for file_name in file_names:
        getH(file_name)







def getValues(keyStr):
    databases = [v for v in glob.glob("{0:}/*.odb".format(wd)) if keyStr in v.split("\\")[-1]]
    for v in databases: 
        extractCurves(v)
    #
    getMax()



# databases = [v for v in glob.glob("{0:}/*.odb".format(wd)) if "aL" in v.split("\\")[-1]]
# for v in databases: 
#     extractCurves(v)


# getMax()


# file_names = databases #glob.glob(os.getcwd()+"/*.odb")


# getHashin()

# # END OF HASHIN DATA EXTRACTION =========================================================================

# exportLoadCurves()








