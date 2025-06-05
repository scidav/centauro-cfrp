# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import optimization
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior



import numpy as np
import os
import glob


session.journalOptions.setValues(replayGeometry=INDEX)



def createTube(nameModel, tubeLength, tubeDiameterOut, tubeThickness,damageLocation):
    if damageLocation=="null":
        dloc = 0.5*tubeLength
    else:
        dloc = damageLocation
    #
    myModel = mdb.Model(name=nameModel, modelType=STANDARD_EXPLICIT)
    if'Model-1' in mdb.models.keys():
        del mdb.models['Model-1']
    s = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, 0.5*(tubeDiameterOut-tubeThickness)))
    p = myModel.Part(name='TUBE', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = myModel.parts['TUBE']
    p.BaseShellExtrude(sketch=s, depth=tubeLength)
    s.unsetPrimaryObject()
    p = myModel.parts['TUBE']
    del myModel.sketches['__profile__']
    # cliCommand("""session.journalOptions.setValues(replayGeometry=INDEX)""")
    leff = tubeDiameterOut+2*max(0.25*tubeDiameterOut,200.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=dloc-0.5*leff)
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=dloc+0.5*leff)
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=dloc)
    p.PartitionFaceByDatumPlane(datumPlane=p.datums[2], faces=p.faces)
    p.PartitionFaceByDatumPlane(datumPlane=p.datums[3], faces=p.faces)
    p.PartitionFaceByDatumPlane(datumPlane=p.datums[4], faces=p.faces)
    p.PartitionFaceByDatumPlane(datumPlane=p.datums[5], faces=p.faces)
    p.PartitionFaceByDatumPlane(datumPlane=p.datums[6], faces=p.faces)
    p = myModel.parts['TUBE']
    p.Set(edges=p.edges[28:29]+p.edges[32:33]+p.edges[34:36], name='Set-BOT')
    p.Set(edges=p.edges[16:17]+p.edges[21:22]+p.edges[25:26]+p.edges[30:31], name='Set-TOP')
    myModel.Material(name='API-5L-GrB')
    myModel.materials['API-5L-GrB'].Elastic(table=((210000.0, 0.3), ))
    myModel.materials['API-5L-GrB'].Plastic(table=((245.0, 0.0), (245.0, 0.18)))
    myModel.HomogeneousShellSection(name='TUBE-Section', 
        preIntegrate=OFF, material='API-5L-GrB', thicknessType=UNIFORM, 
        thickness=tubeThickness, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=5)
    region = p.Set(faces=p.faces, name='TUBE-Body')
    p.SectionAssignment(region=region, 
                        sectionName='TUBE-Section', 
                        offset=0.0, 
                        offsetType=MIDDLE_SURFACE, 
                        offsetField='', 
                        thicknessAssignment=FROM_SECTION)
    p.setMeshControls(regions=p.faces, elemShape=QUAD, technique=FREE, algorithm=ADVANCING_FRONT)#MEDIAL_AXIS, minTransition=ON)# technique=STRUCTURED, minTransition=OFF)
    p.seedPart(size=2.5*tubeThickness, deviationFactor=0.025, minSizeFactor=0.025)
    p.generateMesh()
    a = myModel.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='TUBE-1', part=p, dependent=ON)
    a.ReferencePoint(point=(0.0, np.floor(tubeDiameterOut/110), 0.0))
    a.ReferencePoint(point=(0.0, np.floor(tubeDiameterOut/110), tubeLength))
    refPoints1=(a.referencePoints[4], )
    refPoints2=(a.referencePoints[5], )
    a.Set(referencePoints=refPoints1, name='RP-LOAD')
    a.Set(referencePoints=refPoints2, name='RP-REACTION')
    region4a=a.instances['TUBE-1'].sets['Set-TOP']
    region4b=a.instances['TUBE-1'].sets['Set-BOT']
    region1a=a.sets['RP-REACTION']
    region1b=a.sets['RP-LOAD']
    myModel.RigidBody(name='Constraint-BOT', refPointRegion=region1a, tieRegion=region4a)
    myModel.RigidBody(name='Constraint-TOP', refPointRegion=region1b, tieRegion=region4b)
    a = myModel.rootAssembly
    region = a.sets['RP-LOAD']
    myModel.DisplacementBC(name='BC-LOAD', 
                           createStepName='Initial', region=region, 
                           u1 =SET,    u2=SET,  u3=UNSET, 
                           ur1=UNSET, ur2=SET, ur3=SET, 
                           amplitude=UNSET, distributionType=UNIFORM, 
                           fieldName='', localCsys=None)
    a = myModel.rootAssembly
    region = a.sets['RP-REACTION']
    myModel.DisplacementBC(name='BC-REACTION', 
                           createStepName='Initial', region=region, 
                           u1 =SET,    u2=SET,  u3=SET, 
                           ur1=UNSET, ur2=SET, ur3=SET, 
                           amplitude=UNSET, distributionType=UNIFORM, 
                           fieldName='', localCsys=None)




def setDamage(nameModel, tubeLength, tubeDiameterOut, tubeThickness, damageAngle, damageLocation):
    if damageLocation == "null":
        dloc = 0.5*tubeLength
    else:
        dloc = damageLocation
    #
    myModel = mdb.models[nameModel]
    p = myModel.parts['TUBE']
    t = p.MakeSketchTransform(sketchPlane=p.datums[3], sketchUpEdge=p.edges[35], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, dloc))
    s = myModel.ConstrainedSketch(name='__profile__', sheetSize=4443.78, gridSpacing=111.09, transform=t)
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    damageRadius = 0.5*tubeDiameterOut*np.sin(np.radians(damageAngle))
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, damageRadius))
    p.CutExtrude(sketchPlane=p.datums[3], sketchUpEdge=p.edges[35], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=ON)
    s.unsetPrimaryObject()
    del myModel.sketches['__profile__']
    p.seedEdgeBySize(edges=p.edges[0:1]+p.edges[5:6]+p.edges[16:17]+p.edges[25:26], size=0.3*tubeThickness, deviationFactor=0.01, minSizeFactor=0.01, constraint=FINER)
    p.generateMesh()



def setRepair(nameModel, tubeLength, tubeDiameterOut, tubeThickness, damageAngle, repairThickness, cfrp_dict, damageLocation):
    if damageLocation == "null":
        dloc = 0.5*tubeLength
    else:
        dloc = damageLocation
    #
    myModel = mdb.models[nameModel]
    plyThickness = cfrp_dict['Ply.Thickness']
    s1 = myModel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d1, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, 0.5*(tubeDiameterOut+repairThickness)))
    p = myModel.Part(name='REPAIR', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = myModel.parts['REPAIR']
    damageRadius = 0.5*tubeDiameterOut*np.sin(np.radians(damageAngle))
    leff = 2.0*(damageRadius+max(0.25*tubeDiameterOut, 200.0))
    p.BaseShellExtrude(sketch=s1, depth=leff)
    s1.unsetPrimaryObject()
    p.Surface(side1Faces=p.faces, name='REPAIR-OUT')
    p.Surface(side2Faces=p.faces, name='REPAIR-IN')
    p.Set(faces=p.faces, name='Set-REPAIR')
    matRep = myModel.Material(name=cfrp_dict['Name'])
    me = cfrp_dict['Elastic']
    matRep.Elastic(type=me[0], table=(me[1], ))
    matRep.HashinDamageInitiation(table=(cfrp_dict["Hashin.Initiation"], ))
    h1, h2, h3, h4 = cfrp_dict["Hashin.Stabilization"]
    matRep.hashinDamageInitiation.DamageStabilization(fiberTensileCoeff =h1, fiberCompressiveCoeff =h2, 
                                                      matrixTensileCoeff=h3, matrixCompressiveCoeff=h4)
    matRep.hashinDamageInitiation.DamageEvolution(type=ENERGY, table=((10.9, 0.95, 17.6, 11.55), ))
    layupOrientation = None
    region1=p.sets['Set-REPAIR']
    compositeLayup = myModel.parts['REPAIR'].CompositeLayup(
        name='CompositeLayup-1', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    #===================================================
    # LOOP FOR NUMBER OF PLYS
    #===================================================
    nPlyes = np.ceil(repairThickness/plyThickness)
    for i in range(int(nPlyes)):
        compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(i+1),
            region=region1, material=cfrp_dict['Name'], thicknessType=SPECIFY_THICKNESS,
            thickness=plyThickness, orientationType=SPECIFY_ORIENT, orientationValue=90.0,
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
    # ==================================================
    # END OF LOOP
    # ==================================================
    p.seedPart(size=2*tubeThickness, deviationFactor=0.05, minSizeFactor=0.05)
    p.setMeshControls(regions=p.faces, elemShape=QUAD)
    p.generateMesh()
    p = myModel.parts['TUBE']
    p.Surface(side1Faces=p.faces, name='TUBE-OUT')
    a = myModel.rootAssembly
    a.regenerate()
    p = myModel.parts['REPAIR']
    a.Instance(name='REPAIR-1', part=p, dependent=ON)
    a.translate(instanceList=('REPAIR-1', ), vector=(0.0, 0.0, dloc-0.5*leff))
    region1=a.instances['TUBE-1'].surfaces['TUBE-OUT']
    region2=a.instances['REPAIR-1'].surfaces['REPAIR-IN']
    myModel.Tie(name='Repair', master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, 
        thickness=ON)



def simulation(nameModel, nCpus):
    myModel = mdb.models[nameModel]
    myModel.StaticStep(name='LoadStep', previous='Initial', 
        maxNumInc=10000000, stabilizationMagnitude=0.0002, 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION, 
        continueDampingFactors=False, adaptiveDampingRatio=0.05, 
        initialInc=0.01, minInc=1e-05, maxInc=0.01, nlgeom=ON)
    myModel.boundaryConditions['BC-LOAD'].setValuesInStep(
        stepName='LoadStep', u3=10.0)
    myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 
        'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
        'HSNFTCRT', 'HSNFCCRT', 'HSNMTCRT', 'HSNMCCRT', 'COORD'))
    mdb.Job(name='J-'+nameModel, model=nameModel, description='', type=ANALYSIS, atTime=None, 
        waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=nCpus, 
        numDomains=nCpus, numGPUs=1)
    mdb.jobs['J-'+nameModel].writeInput(consistencyChecking=OFF)



def correctDamagedModel(nameModel):
    myModel = mdb.models[nameModel]
    myModel.steps['LoadStep'].setValues(stabilizationMethod=NONE, 
                                        continueDampingFactors=False, 
                                        adaptiveDampingRatio=None, 
                                        initialInc=0.01, 
                                        minInc=3.0e-04, 
                                        maxInc=0.03)
    mdb.jobs['J-'+nameModel].writeInput(consistencyChecking=OFF)
    


def correctRepairedModel(nameModel):
    myModel = mdb.models[nameModel]
    a = myModel.rootAssembly
    p = myModel.parts['TUBE']
    p.Surface(side1Faces=p.faces[0:5]+p.faces[7:8]+p.faces[12:14], name='TUBE-OUT')
    p.generateMesh()
    a.regenerate()
    a = myModel.rootAssembly
    myModel.steps['LoadStep'].setValues(stabilizationMethod=NONE, 
                                        continueDampingFactors=False, 
                                        adaptiveDampingRatio=None, 
                                        initialInc=0.01, 
                                        minInc=3.0e-04, 
                                        maxInc=0.03)
    mdb.jobs['J-'+nameModel].writeInput(consistencyChecking=OFF)



def writeBatchFile(nCpus):
    with open(os.getcwd()+"/01 - submit_files.bat","w+") as batchFile:
        batchFile.write('pushd %~dp0\n')
        #
        for jobName in mdb.jobs.keys():
            batchFile.write('call abq2020 interactive job="{0:}" cpus={1:}\n'.format(jobName,nCpus))
        #
        batchFile.write("popd\n")
    #
    batchFile.close()


def saveCAE(nameCaeOut):
    mdb.saveAs(pathName='{0:}/{1:}'.format(os.getcwd(), nameCaeOut))


# ================================================================= END OF MODELING SNIPPETS =================================================================

# ============================================================== START OF POSPROCESSING SNIPPETS =============================================================

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

# END OF DATA WRITTING FUNCTIONS ============================================


# BEGIN OF HASHIN DATA EXTRACTION ===========================================

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


