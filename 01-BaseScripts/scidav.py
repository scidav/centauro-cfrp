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




# =======================================================================================
#  END OF LAYERED PROGRAMMING
# =======================================================================================


composite_MA = { 
                    'Name'                  : 'COMPOSITE-Macseal',
                    'Elastic'               : [LAMINA,[60543.0, 27080.0, 0.078, 2601.0, 1571.48, 1571.48]],
                    'Hashin.Initiation'     : [965.560, 345.96, 212.309, 133.590, 41.621, 38.819],
                    #'Hashin.Evolution'      : [ENERGY,[10.9, 0.95, 17.6, 11.55]],
                    'Hashin.Stabilization'  : [0.001, 0.001, 0.001, 0.001],
                    'Ply.Thickness' : 0.54 }

composite_TK = { 
            'Name'                  : 'COMPOSITE-Tecnofink',
            'Elastic'               : [LAMINA,[52.390E3, 5.259E3, 0.316, 2068.0, 1864.1, 1864.1]],
            'Hashin.Initiation'     : [598.407, 212.484, 18.609, 67.401, 42.816, 19.160],
            #'Hashin.Evolution'      : [ENERGY,[10.9, 0.95, 17.6, 11.55]],
            'Hashin.Stabilization'  : [0.001, 0.001, 0.001, 0.001],
            'Ply.Thickness' : 0.9 }

composite_IT = { 
            'Name'                  : 'COMPOSITE-HITA',
            'Elastic'               : [LAMINA,[24584.0, 16378.0, 0.322, 4686.0, 1124.8, 1124.8]],
            'Hashin.Initiation'     : [308.530, 174.450, 193.651, 133.748, 68.112, 76.777],
            #'Hashin.Evolution'      : [ENERGY,[10.9, 0.95, 17.6, 11.55]],
            'Hashin.Stabilization'  : [0.001, 0.001, 0.001, 0.001],
            'Ply.Thickness': 1.0 }



# Here declare the material used for repair if repaired models are being simulated
cfrp_mat = composite_MA


# doe for damage models with aL damage location
doe = [['M01-D15-aL12', 3200.00, 246.00, 12.3000, 15.00, 1600.00],
        ['M01-D30-aL12', 3200.00, 246.00, 12.3000, 30.00, 1600.00],
        ['M01-D45-aL12', 3200.00, 246.00, 12.3000, 45.00, 1600.00],
        ['M01-D60-aL12', 3200.00, 246.00, 12.3000, 60.00, 1600.00],
        ['M02-D15-aL12', 3473.00, 178.00, 9.8889, 15.00, 1736.50],
        ['M02-D30-aL12', 3473.00, 178.00, 9.8889, 30.00, 1736.50],
        ['M02-D45-aL12', 3473.00, 178.00, 9.8889, 45.00, 1736.50],
        ['M02-D60-aL12', 3473.00, 178.00, 9.8889, 60.00, 1736.50],
        ['M03-D15-aL12', 4397.00, 169.00, 4.8286, 15.00, 2198.50],
        ['M03-D30-aL12', 4397.00, 169.00, 4.8286, 30.00, 2198.50],
        ['M03-D45-aL12', 4397.00, 169.00, 4.8286, 45.00, 2198.50],
        ['M03-D60-aL12', 4397.00, 169.00, 4.8286, 60.00, 2198.50],
        ['M04-D15-aL12', 9788.00, 301.00, 14.3333, 15.00, 4894.00],
        ['M04-D30-aL12', 9788.00, 301.00, 14.3333, 30.00, 4894.00],
        ['M04-D45-aL12', 9788.00, 301.00, 14.3333, 45.00, 4894.00],
        ['M04-D60-aL12', 9788.00, 301.00, 14.3333, 60.00, 4894.00],
        ['M05-D15-aL12', 7414.00, 190.00, 6.5517, 15.00, 3707.00],
        ['M05-D30-aL12', 7414.00, 190.00, 6.5517, 30.00, 3707.00],
        ['M05-D45-aL12', 7414.00, 190.00, 6.5517, 45.00, 3707.00],
        ['M05-D60-aL12', 7414.00, 190.00, 6.5517, 60.00, 3707.00],
        ['M06-D15-aL12', 7057.00, 155.00, 4.4286, 15.00, 3528.50],
        ['M06-D30-aL12', 7057.00, 155.00, 4.4286, 30.00, 3528.50],
        ['M06-D45-aL12', 7057.00, 155.00, 4.4286, 45.00, 3528.50],
        ['M06-D60-aL12', 7057.00, 155.00, 4.4286, 60.00, 3528.50],
        ['M07-D15-aL12', 5827.00, 112.00, 4.4800, 15.00, 2913.50],
        ['M07-D30-aL12', 5827.00, 112.00, 4.4800, 30.00, 2913.50],
        ['M07-D45-aL12', 5827.00, 112.00, 4.4800, 45.00, 2913.50],
        ['M07-D60-aL12', 5827.00, 112.00, 4.4800, 60.00, 2913.50],
        ['M01-D15-aL13', 3200.00, 246.00, 12.3000, 15.00, 1066.67],
        ['M01-D30-aL13', 3200.00, 246.00, 12.3000, 30.00, 1066.67],
        ['M01-D45-aL13', 3200.00, 246.00, 12.3000, 45.00, 1066.67],
        ['M01-D60-aL13', 3200.00, 246.00, 12.3000, 60.00, 1066.67],
        ['M02-D15-aL13', 3473.00, 178.00, 9.8889, 15.00, 1157.67],
        ['M02-D30-aL13', 3473.00, 178.00, 9.8889, 30.00, 1157.67],
        ['M02-D45-aL13', 3473.00, 178.00, 9.8889, 45.00, 1157.67],
        ['M02-D60-aL13', 3473.00, 178.00, 9.8889, 60.00, 1157.67],
        ['M03-D15-aL13', 4397.00, 169.00, 4.8286, 15.00, 1465.67],
        ['M03-D30-aL13', 4397.00, 169.00, 4.8286, 30.00, 1465.67],
        ['M03-D45-aL13', 4397.00, 169.00, 4.8286, 45.00, 1465.67],
        ['M03-D60-aL13', 4397.00, 169.00, 4.8286, 60.00, 1465.67],
        ['M04-D15-aL13', 9788.00, 301.00, 14.3333, 15.00, 3262.67],
        ['M04-D30-aL13', 9788.00, 301.00, 14.3333, 30.00, 3262.67],
        ['M04-D45-aL13', 9788.00, 301.00, 14.3333, 45.00, 3262.67],
        ['M04-D60-aL13', 9788.00, 301.00, 14.3333, 60.00, 3262.67],
        ['M05-D15-aL13', 7414.00, 190.00, 6.5517, 15.00, 2471.33],
        ['M05-D30-aL13', 7414.00, 190.00, 6.5517, 30.00, 2471.33],
        ['M05-D45-aL13', 7414.00, 190.00, 6.5517, 45.00, 2471.33],
        ['M05-D60-aL13', 7414.00, 190.00, 6.5517, 60.00, 2471.33],
        ['M06-D15-aL13', 7057.00, 155.00, 4.4286, 15.00, 2352.33],
        ['M06-D30-aL13', 7057.00, 155.00, 4.4286, 30.00, 2352.33],
        ['M06-D45-aL13', 7057.00, 155.00, 4.4286, 45.00, 2352.33],
        ['M06-D60-aL13', 7057.00, 155.00, 4.4286, 60.00, 2352.33],
        ['M07-D15-aL13', 5827.00, 112.00, 4.4800, 15.00, 1942.33],
        ['M07-D30-aL13', 5827.00, 112.00, 4.4800, 30.00, 1942.33],
        ['M07-D45-aL13', 5827.00, 112.00, 4.4800, 45.00, 1942.33],
        ['M07-D60-aL13', 5827.00, 112.00, 4.4800, 60.00, 1942.33],
        ['M01-D15-aL16', 3200.00, 246.00, 12.3000, 15.00, 533.33],
        ['M01-D30-aL16', 3200.00, 246.00, 12.3000, 30.00, 533.33],
        ['M01-D45-aL16', 3200.00, 246.00, 12.3000, 45.00, 533.33],
        ['M01-D60-aL16', 3200.00, 246.00, 12.3000, 60.00, 533.33],
        ['M02-D15-aL16', 3473.00, 178.00, 9.8889, 15.00, 578.83],
        ['M02-D30-aL16', 3473.00, 178.00, 9.8889, 30.00, 578.83],
        ['M02-D45-aL16', 3473.00, 178.00, 9.8889, 45.00, 578.83],
        ['M02-D60-aL16', 3473.00, 178.00, 9.8889, 60.00, 578.83],
        ['M03-D15-aL16', 4397.00, 169.00, 4.8286, 15.00, 732.83],
        ['M03-D30-aL16', 4397.00, 169.00, 4.8286, 30.00, 732.83],
        ['M03-D45-aL16', 4397.00, 169.00, 4.8286, 45.00, 732.83],
        ['M03-D60-aL16', 4397.00, 169.00, 4.8286, 60.00, 732.83],
        ['M04-D15-aL16', 9788.00, 301.00, 14.3333, 15.00, 1631.33],
        ['M04-D30-aL16', 9788.00, 301.00, 14.3333, 30.00, 1631.33],
        ['M04-D45-aL16', 9788.00, 301.00, 14.3333, 45.00, 1631.33],
        ['M04-D60-aL16', 9788.00, 301.00, 14.3333, 60.00, 1631.33],
        ['M05-D15-aL16', 7414.00, 190.00, 6.5517, 15.00, 1235.67],
        ['M05-D30-aL16', 7414.00, 190.00, 6.5517, 30.00, 1235.67],
        ['M05-D45-aL16', 7414.00, 190.00, 6.5517, 45.00, 1235.67],
        ['M05-D60-aL16', 7414.00, 190.00, 6.5517, 60.00, 1235.67],
        ['M06-D15-aL16', 7057.00, 155.00, 4.4286, 15.00, 1176.17],
        ['M06-D30-aL16', 7057.00, 155.00, 4.4286, 30.00, 1176.17],
        ['M06-D45-aL16', 7057.00, 155.00, 4.4286, 45.00, 1176.17],
        ['M06-D60-aL16', 7057.00, 155.00, 4.4286, 60.00, 1176.17],
        ['M07-D15-aL16', 5827.00, 112.00, 4.4800, 15.00, 971.17],
        ['M07-D30-aL16', 5827.00, 112.00, 4.4800, 30.00, 971.17],
        ['M07-D45-aL16', 5827.00, 112.00, 4.4800, 45.00, 971.17],
        ['M07-D60-aL16', 5827.00, 112.00, 4.4800, 60.00, 971.17]]




for m in doe:
    # nameModel1, tubeLength1, tubeDiameterOut1, tubeThickness1, damageAngle1, damageLocation1, repairThickness1 = m[0], m[1], m[2], m[3], m[4], m[5], m[6]#  for repair with wi/aL
    nameModel1, tubeLength1, tubeDiameterOut1, tubeThickness1, damageAngle1, damageLocation1  = m[0], m[1], m[2], m[3], m[4], m[5] #  for damage wi/aL
    createTube(nameModel1, tubeLength1, tubeDiameterOut1, tubeThickness1,damageLocation1)
    setDamage(nameModel1, tubeLength1, tubeDiameterOut1, tubeThickness1, damageAngle1, damageLocation1)
    # setRepair(nameModel1, tubeLength1, tubeDiameterOut1, tubeThickness1, damageAngle1, repairThickness1, cfrp_mat, damageLocation1) # uncomment if repair
    simulation(nameModel1, 20)
    # correctRepairedModel(nameModel1) # uncomment if repair
    correctDamagedModel(nameModel1)


# ======================================================================================================================================================
# Rememeber that for optimization, the Leff was defined as the ceil integer of the estimated Leff parameter and the one near to 50 to 100 multiple
# ======================================================================================================================================================


writeBatchFile(20) # set rthe nCPUs
saveCAE("Damage-AxialLocationEffect-202506-V0") #set a name for the cae file




