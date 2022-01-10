#!/usr/bin/env python

# load up the opencmiss library of functions into an object called iron.
from opencmiss.iron import iron

import re
import numpy as np
import math
import random
import os
from inputs import *
from modules import select
from forcebased_k_new import alpha_dis, k5r_dis_cont, myosin_dis

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# -------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
# -------------------------------------------------------------------------------------------------

# Two Dimensional Coordinate System
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(CoordinateSystemUserNumber)
coordinateSystem.dimension = 2
coordinateSystem.CreateFinish()

# -------------------------------------------------------------------------------------------------
# REGION
# -------------------------------------------------------------------------------------------------

# Start Region
region = iron.Region()
region.CreateStart(RegionUserNumber, iron.WorldRegion)
region.label = "membrane"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# -------------------------------------------------------------------------------------------------
# BASIS
# -------------------------------------------------------------------------------------------------

# Simplex Basis
basis = iron.Basis()
basis.CreateStart(BasisUserNumber)
basis.type = iron.BasisTypes.SIMPLEX
basis.numberOfXi = numberOfXi
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * numberOfXi
# basis.quadratureNumberOfGaussXi = [2]*2
basis.CreateFinish()

# -------------------------------------------------------------------------------------------------
# MESH
# -------------------------------------------------------------------------------------------------

    # ___________________________________________________________________________________
    #                              Initializing Node File

# Inputing node file
node_file = open(nodeinput, 'r')

# Reading the values of the first line
num_nodes, num_coords, num_attributes, boundary_markers = node_file.readline().split()

# Converting the numbers to integers
num_nodes = int(num_nodes)
num_coords = int(num_coords)
num_attributes = int(num_attributes)
boundary_marker = int(boundary_markers)

# Creating variables to store node number & boundary marker
NodeNums = [[0 for m in range(2)] for n in range(num_nodes)]

# Creating variable to store x, y and z coordinates
NodeCoords = [[0 for m in range(num_coords)] for n in range(num_nodes)]

# Reading data from nodefiler
for i in range(num_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeNums[i][1] = node_file.readline().split()
    # Converting from string to int/float
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeNums[i][1] = int(NodeNums[i][1])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])
    
# The number of Boundary Nodes
Num_BoundaryNodes = [row[1] for row in NodeNums].count(1)

# Closing the file
node_file.close()

# ___________________________________________________________________________________
#                              Initializing Element File

# Inputing element file
eleminput = eleminput.replace("\r","")
elem_file = open(eleminput, 'r')

# Reading the values of the first line
num_elements, nodes_per_ele, ele_attributes = elem_file.readline().split()

# Converting the numbers to integers
num_elements = int(num_elements)
nodes_per_ele = int(nodes_per_ele)
ele_attributes = int(ele_attributes)

# Creating variable to store element map
ElemMap = [[0 for x in range(3)] for y in range(num_elements)]
Elemindex = [[0 for m in range(1)] for n in range(num_elements)]
# Reading data from elemfile
for i in range(num_elements):

    Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2] = elem_file.readline().split()

    # Converting from string to int
    Elemindex[i][0] = int(Elemindex[i][0])
    ElemMap[i][0] = int(ElemMap[i][0])
    ElemMap[i][1] = int(ElemMap[i][1])
    ElemMap[i][2] = int(ElemMap[i][2])

# Closing the file
elem_file.close()

###############################################################################
# Incenter node and element file
###############################################################################

node_incenter = np.loadtxt(nodeincenter)
ele_incenter = np.loadtxt(elemincenter)
num_nodes_incenter = int(node_incenter[0,0])
# Initialise Nodes
nodes = iron.Nodes()
nodes.CreateStart(region, num_nodes_incenter)
nodes.CreateFinish()

# Initialise Mesh
mesh = iron.Mesh()
mesh.CreateStart(MeshUserNumber, region, int(node_incenter[0,1]))
mesh.NumberOfElementsSet(int(ele_incenter[0,0]))
mesh.NumberOfComponentsSet(1)

# Initialise Elements
meshElements = iron.MeshElements()
meshElements.CreateStart(mesh, 1, basis)
for i in range(int(ele_incenter[0,0])):

    meshElements.NodesSet(int(ele_incenter[i+1,0]),\
                         [int(ele_incenter[i+1,1]),\
                          int(ele_incenter[i+1,2]),\
                          int(ele_incenter[i+1,3])])
meshElements.CreateFinish()

# Finilise Mesh
mesh.CreateFinish()

# -------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
# -------------------------------------------------------------------------------------------------

# Parallelization
decomposition = iron.Decomposition()
decomposition.CreateStart(DecompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# -------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
# -------------------------------------------------------------------------------------------------

# Geometric Field
geometricField = iron.Field()
geometricField.CreateStart(GeometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)

geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
geometricField.CreateFinish()

# Update Geometric Field from customized mesh
for i in range(num_nodes_incenter):
    i=i+1
    node = int(node_incenter[i,0])
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain == computationalNodeNumber:
        nodex = node_incenter[i,1]
        nodey = node_incenter[i,2]
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 1, nodex)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                iron.FieldParameterSetTypes.VALUES,
                                                1, 1, node, 2, nodey)


# Update Geometric Field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                       iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)


#_________________________________________________________________________________________________
# REACTION DIFFUSION
#_________________________________________________________________________________________________

#-------------------------------------------------------------------------------------------------
# E_trans EQUATIONS (THE FIRST PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - E_trans
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
E_trans_EquationsSetField         = iron.Field()
E_trans_EquationsSet              = iron.EquationsSet()
E_trans_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
E_trans_EquationsSet.CreateStart(E_trans_EquationsSetUserNumber,region,
                                   geometricField,E_trans_EquationsSetSpecification,
                                   E_trans_EquationsSetFieldUserNumber,E_trans_EquationsSetField)
E_trans_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - E_trans
     #--------------------------------------------------------------------------------------------

#Start and label E_trans Field Dependent Field
E_trans_Field = iron.Field()
E_trans_EquationsSet.DependentCreateStart(E_trans_FieldUserNumber, E_trans_Field)
E_trans_Field.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans Field")
E_trans_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"E_trans Field DELUDELN")

#Finish E_trans Field Dependent Field
E_trans_EquationsSet.DependentCreateFinish()

#Initialise E_trans Field Dependent field
E_trans_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_E_trans)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - E_trans
     #--------------------------------------------------------------------------------------------

E_trans_MaterialsField = iron.Field()
E_trans_EquationsSet.MaterialsCreateStart(E_trans_MaterialsFieldUserNumber,E_trans_MaterialsField)

E_trans_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans Materials Field")
E_trans_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
E_trans_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_E_trans)

#Diffusion Coefficient in Y
E_trans_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_E_trans)

#Storage Coefficient
E_trans_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - E_trans
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
E_trans_SourceField = iron.Field()
E_trans_EquationsSet.SourceCreateStart(E_trans_SourceFieldUserNumber, E_trans_SourceField)
E_trans_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans Source Field")
E_trans_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
E_trans_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update E_trans Source Field
E_trans_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
E_trans_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update E_trans Dependent Field
E_trans_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
E_trans_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)




#-------------------------------------------------------------------------------------------------
# actin_1 EQUATIONS (THE SECOND PROTEIN)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - actin_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
actin_1_EquationsSetField         = iron.Field()
actin_1_EquationsSet              = iron.EquationsSet()
actin_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
actin_1_EquationsSet.CreateStart(actin_1_EquationsSetUserNumber,region,
                                   geometricField,actin_1_EquationsSetSpecification,
                                   actin_1_EquationsSetFieldUserNumber,actin_1_EquationsSetField)
actin_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - actin_1
     #--------------------------------------------------------------------------------------------

#Start and label actin_1 Field Dependent Field
actin_1_Field = iron.Field()
actin_1_EquationsSet.DependentCreateStart(actin_1_FieldUserNumber, actin_1_Field)
actin_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"actin_1 Field")
actin_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"actin_1 Field DELUDELN")

#Finish actin_1 Field Dependent Field
actin_1_EquationsSet.DependentCreateFinish()

#Initialise actin_1 Field Dependent field
actin_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_actin_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - actin_1
     #--------------------------------------------------------------------------------------------

actin_1_MaterialsField = iron.Field()
actin_1_EquationsSet.MaterialsCreateStart(actin_1_MaterialsFieldUserNumber,actin_1_MaterialsField)

actin_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_1 Materials Field")
actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
actin_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_actin)

#Diffusion Coefficient in Y
actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_actin)

#Storage Coefficient
actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - actin_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
actin_1_SourceField = iron.Field()
actin_1_EquationsSet.SourceCreateStart(actin_1_SourceFieldUserNumber, actin_1_SourceField)
actin_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_1 Source Field")
actin_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
actin_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update actin_1 Source Field
actin_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
actin_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update actin_1 Dependent Field
actin_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
actin_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# actin_a_1 EQUATIONS (Activated actin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - actin_a_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
actin_a_1_EquationsSetField         = iron.Field()
actin_a_1_EquationsSet              = iron.EquationsSet()
actin_a_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
actin_a_1_EquationsSet.CreateStart(actin_a_1_EquationsSetUserNumber,region,
                                   geometricField,actin_a_1_EquationsSetSpecification,
                                   actin_a_1_EquationsSetFieldUserNumber,actin_a_1_EquationsSetField)
actin_a_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - actin_a_1
     #--------------------------------------------------------------------------------------------

#Start and label actin_a_1 Field Dependent Field
actin_a_1_Field = iron.Field()
actin_a_1_EquationsSet.DependentCreateStart(actin_a_1_FieldUserNumber, actin_a_1_Field)
actin_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"actin_a_1 Field")
actin_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"actin_a_1 Field DELUDELN")

#Finish actin_a_1 Field Dependent Field
actin_a_1_EquationsSet.DependentCreateFinish()

#Initialise actin_a_1 Field Dependent field
actin_a_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_actin_a_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - actin_a_1
     #--------------------------------------------------------------------------------------------

actin_a_1_MaterialsField = iron.Field()
actin_a_1_EquationsSet.MaterialsCreateStart(actin_a_1_MaterialsFieldUserNumber,actin_a_1_MaterialsField)

actin_a_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_a_1 Materials Field")
actin_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
actin_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
actin_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
actin_a_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
actin_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_actin_a)

#Diffusion Coefficient in Y
actin_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_actin_a)

#Storage Coefficient
actin_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - actin_a_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
actin_a_1_SourceField = iron.Field()
actin_a_1_EquationsSet.SourceCreateStart(actin_a_1_SourceFieldUserNumber, actin_a_1_SourceField)
actin_a_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_a_1 Source Field")
actin_a_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
actin_a_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update actin_a_1 Source Field
actin_a_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
actin_a_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update actin_a_1 Dependent Field
actin_a_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
actin_a_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# Rac_1 EQUATIONS
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - Rac_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
Rac_1_EquationsSetField         = iron.Field()
Rac_1_EquationsSet              = iron.EquationsSet()
Rac_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
Rac_1_EquationsSet.CreateStart(Rac_1_EquationsSetUserNumber,region,
                                   geometricField,Rac_1_EquationsSetSpecification,
                                   Rac_1_EquationsSetFieldUserNumber,Rac_1_EquationsSetField)
Rac_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - Rac_1
     #--------------------------------------------------------------------------------------------

#Start and label Rac_1 Field Dependent Field
Rac_1_Field = iron.Field()
Rac_1_EquationsSet.DependentCreateStart(Rac_1_FieldUserNumber, Rac_1_Field)
Rac_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_1 Field")
Rac_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"Rac_1 Field DELUDELN")

#Finish Rac_1 Field Dependent Field
Rac_1_EquationsSet.DependentCreateFinish()

#Initialise Rac_1 Field Dependent field
Rac_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_Rac_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - Rac_1
     #--------------------------------------------------------------------------------------------

Rac_1_MaterialsField = iron.Field()
Rac_1_EquationsSet.MaterialsCreateStart(Rac_1_MaterialsFieldUserNumber,Rac_1_MaterialsField)

Rac_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_1 Materials Field")
Rac_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
Rac_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
Rac_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
Rac_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
Rac_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_Rac)

#Diffusion Coefficient in Y
Rac_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_Rac)

#Storage Coefficient
Rac_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)

     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - Rac_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
Rac_1_SourceField = iron.Field()
Rac_1_EquationsSet.SourceCreateStart(Rac_1_SourceFieldUserNumber, Rac_1_SourceField)
Rac_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_1 Source Field")
Rac_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
Rac_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update Rac_1 Source Field
Rac_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
Rac_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update Rac_1 Dependent Field
Rac_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
Rac_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# Rac_a_1 EQUATIONS (Activated Rac)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - Rac_a_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
Rac_a_1_EquationsSetField         = iron.Field()
Rac_a_1_EquationsSet              = iron.EquationsSet()
Rac_a_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
Rac_a_1_EquationsSet.CreateStart(Rac_a_1_EquationsSetUserNumber,region,
                                   geometricField,Rac_a_1_EquationsSetSpecification,
                                   Rac_a_1_EquationsSetFieldUserNumber,Rac_a_1_EquationsSetField)
Rac_a_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - Rac_a_1
     #--------------------------------------------------------------------------------------------

#Start and label Rac_a_1 Field Dependent Field
Rac_a_1_Field = iron.Field()
Rac_a_1_EquationsSet.DependentCreateStart(Rac_a_1_FieldUserNumber, Rac_a_1_Field)
Rac_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_a_1 Field")
Rac_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"Rac_a_1 Field DELUDELN")

#Finish Rac_a_1 Field Dependent Field
Rac_a_1_EquationsSet.DependentCreateFinish()

#Initialise Rac_a_1 Field Dependent field
Rac_a_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_Rac_a_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - Rac_a_1
     #--------------------------------------------------------------------------------------------

Rac_a_1_MaterialsField = iron.Field()
Rac_a_1_EquationsSet.MaterialsCreateStart(Rac_a_1_MaterialsFieldUserNumber,Rac_a_1_MaterialsField)

Rac_a_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_a_1 Materials Field")
Rac_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
Rac_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
Rac_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
Rac_a_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
Rac_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_Rac_a)

#Diffusion Coefficient in Y
Rac_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_Rac_a)

#Storage Coefficient
Rac_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - Rac_a_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
Rac_a_1_SourceField = iron.Field()
Rac_a_1_EquationsSet.SourceCreateStart(Rac_a_1_SourceFieldUserNumber, Rac_a_1_SourceField)
Rac_a_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_a_1 Source Field")
Rac_a_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
Rac_a_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update Rac_a_1 Source Field
Rac_a_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
Rac_a_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update Rac_a_1 Dependent Field
Rac_a_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
Rac_a_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# RhoA_1 EQUATIONS (RhoA_1)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - RhoA_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
RhoA_1_EquationsSetField         = iron.Field()
RhoA_1_EquationsSet              = iron.EquationsSet()
RhoA_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
RhoA_1_EquationsSet.CreateStart(RhoA_1_EquationsSetUserNumber,region,
                                   geometricField,RhoA_1_EquationsSetSpecification,
                                   RhoA_1_EquationsSetFieldUserNumber,RhoA_1_EquationsSetField)
RhoA_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - RhoA_1
     #--------------------------------------------------------------------------------------------

#Start and label RhoA_1 Field Dependent Field
RhoA_1_Field = iron.Field()
RhoA_1_EquationsSet.DependentCreateStart(RhoA_1_FieldUserNumber, RhoA_1_Field)
RhoA_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_1 Field")
RhoA_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"RhoA_1 Field DELUDELN")

#Finish RhoA_1 Field Dependent Field
RhoA_1_EquationsSet.DependentCreateFinish()

#Initialise RhoA_1 Field Dependent field
RhoA_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_RhoA_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - RhoA_1
     #--------------------------------------------------------------------------------------------

RhoA_1_MaterialsField = iron.Field()
RhoA_1_EquationsSet.MaterialsCreateStart(RhoA_1_MaterialsFieldUserNumber,RhoA_1_MaterialsField)

RhoA_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_1 Materials Field")
RhoA_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
RhoA_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_RhoA)

#Diffusion Coefficient in Y
RhoA_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_RhoA)

#Storage Coefficient
RhoA_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - RhoA_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
RhoA_1_SourceField = iron.Field()
RhoA_1_EquationsSet.SourceCreateStart(RhoA_1_SourceFieldUserNumber, RhoA_1_SourceField)
RhoA_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_1 Source Field")
RhoA_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
RhoA_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update RhoA_1 Source Field
RhoA_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
RhoA_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update RhoA_1 Dependent Field
RhoA_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
RhoA_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# RhoA_a_1 EQUATIONS (Activated RhoA)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - RhoA_a_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
RhoA_a_1_EquationsSetField         = iron.Field()
RhoA_a_1_EquationsSet              = iron.EquationsSet()
RhoA_a_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
RhoA_a_1_EquationsSet.CreateStart(RhoA_a_1_EquationsSetUserNumber,region,
                                   geometricField,RhoA_a_1_EquationsSetSpecification,
                                   RhoA_a_1_EquationsSetFieldUserNumber,RhoA_a_1_EquationsSetField)
RhoA_a_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - RhoA_a_1
     #--------------------------------------------------------------------------------------------

#Start and label RhoA_a_1 Field Dependent Field
RhoA_a_1_Field = iron.Field()
RhoA_a_1_EquationsSet.DependentCreateStart(RhoA_a_1_FieldUserNumber, RhoA_a_1_Field)
RhoA_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_a_1 Field")
RhoA_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"RhoA_a_1 Field DELUDELN")

#Finish RhoA_a_1 Field Dependent Field
RhoA_a_1_EquationsSet.DependentCreateFinish()

#Initialise RhoA_a_1 Field Dependent field
RhoA_a_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_RhoA_a_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - RhoA_a_1
     #--------------------------------------------------------------------------------------------

RhoA_a_1_MaterialsField = iron.Field()
RhoA_a_1_EquationsSet.MaterialsCreateStart(RhoA_a_1_MaterialsFieldUserNumber,RhoA_a_1_MaterialsField)

RhoA_a_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_a_1 Materials Field")
RhoA_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_a_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
RhoA_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_RhoA_a)

#Diffusion Coefficient in Y
RhoA_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_RhoA_a)

#Storage Coefficient
RhoA_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - RhoA_a_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
RhoA_a_1_SourceField = iron.Field()
RhoA_a_1_EquationsSet.SourceCreateStart(RhoA_a_1_SourceFieldUserNumber, RhoA_a_1_SourceField)
RhoA_a_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_a_1 Source Field")
RhoA_a_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
RhoA_a_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update RhoA_a_1 Source Field
RhoA_a_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
RhoA_a_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update RhoA_a_1 Dependent Field
RhoA_a_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
RhoA_a_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# myo_1 EQUATIONS (Myosin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - myo_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
myo_1_EquationsSetField         = iron.Field()
myo_1_EquationsSet              = iron.EquationsSet()
myo_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
myo_1_EquationsSet.CreateStart(myo_1_EquationsSetUserNumber,region,
                                   geometricField,myo_1_EquationsSetSpecification,
                                   myo_1_EquationsSetFieldUserNumber,myo_1_EquationsSetField)
myo_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - myo_1
     #--------------------------------------------------------------------------------------------

#Start and label myo_1 Field Dependent Field
myo_1_Field = iron.Field()
myo_1_EquationsSet.DependentCreateStart(myo_1_FieldUserNumber, myo_1_Field)
myo_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"myo_1 Field")
myo_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"myo_1 Field DELUDELN")

#Finish myo_1 Field Dependent Field
myo_1_EquationsSet.DependentCreateFinish()

#Initialise myo_1 Field Dependent field
myo_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_myo_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - myo_1
     #--------------------------------------------------------------------------------------------

myo_1_MaterialsField = iron.Field()
myo_1_EquationsSet.MaterialsCreateStart(myo_1_MaterialsFieldUserNumber,myo_1_MaterialsField)

myo_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_1 Materials Field")
myo_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
myo_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
myo_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
myo_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
myo_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_myo)

#Diffusion Coefficient in Y
myo_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_myo)

#Storage Coefficient
myo_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - myo_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
myo_1_SourceField = iron.Field()
myo_1_EquationsSet.SourceCreateStart(myo_1_SourceFieldUserNumber, myo_1_SourceField)
myo_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_1 Source Field")
myo_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
myo_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update myo_1 Source Field
myo_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
myo_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update myo_1 Dependent Field
myo_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
myo_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# myo_a_1 EQUATIONS (Activated Myosin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - myo_a_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
myo_a_1_EquationsSetField         = iron.Field()
myo_a_1_EquationsSet              = iron.EquationsSet()
myo_a_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
myo_a_1_EquationsSet.CreateStart(myo_a_1_EquationsSetUserNumber,region,
                                   geometricField,myo_a_1_EquationsSetSpecification,
                                   myo_a_1_EquationsSetFieldUserNumber,myo_a_1_EquationsSetField)
myo_a_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - myo_a_1
     #--------------------------------------------------------------------------------------------

#Start and label myo_a_1 Field Dependent Field
myo_a_1_Field = iron.Field()
myo_a_1_EquationsSet.DependentCreateStart(myo_a_1_FieldUserNumber, myo_a_1_Field)
myo_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"myo_a_1 Field")
myo_a_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"myo_a_1 Field DELUDELN")

#Finish myo_a_1 Field Dependent Field
myo_a_1_EquationsSet.DependentCreateFinish()

#Initialise myo_a_1 Field Dependent field
myo_a_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_myo_a_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - myo_a_1
     #--------------------------------------------------------------------------------------------

myo_a_1_MaterialsField = iron.Field()
myo_a_1_EquationsSet.MaterialsCreateStart(myo_a_1_MaterialsFieldUserNumber,myo_a_1_MaterialsField)

myo_a_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_a_1 Materials Field")
myo_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
myo_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
myo_a_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
myo_a_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
myo_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_myo_a)

#Diffusion Coefficient in Y
myo_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_myo_a)

#Storage Coefficient
myo_a_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - myo_a_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
myo_a_1_SourceField = iron.Field()
myo_a_1_EquationsSet.SourceCreateStart(myo_a_1_SourceFieldUserNumber, myo_a_1_SourceField)
myo_a_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_a_1 Source Field")
myo_a_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
myo_a_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update myo_a_1 Source Field
myo_a_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
myo_a_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update myo_a_1 Dependent Field
myo_a_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
myo_a_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
#  E_trans_actin_1 EQUATIONS (Bond between E-cad and activated actin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - E_trans_actin_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
E_trans_actin_1_EquationsSetField         = iron.Field()
E_trans_actin_1_EquationsSet              = iron.EquationsSet()
E_trans_actin_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
E_trans_actin_1_EquationsSet.CreateStart(E_trans_actin_1_EquationsSetUserNumber,region,
                                   geometricField,E_trans_actin_1_EquationsSetSpecification,
                                   E_trans_actin_1_EquationsSetFieldUserNumber,E_trans_actin_1_EquationsSetField)
E_trans_actin_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - E_trans_actin_1
     #--------------------------------------------------------------------------------------------

#Start and label E_trans_actin_1 Field Dependent Field
E_trans_actin_1_Field = iron.Field()
E_trans_actin_1_EquationsSet.DependentCreateStart(E_trans_actin_1_FieldUserNumber, E_trans_actin_1_Field)
E_trans_actin_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans_actin_1 Field")
E_trans_actin_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"E_trans_actin_1 Field DELUDELN")

#Finish E_trans_actin_1 Field Dependent Field
E_trans_actin_1_EquationsSet.DependentCreateFinish()

#Initialise E_trans_actin_1 Field Dependent field
E_trans_actin_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_E_trans_actin_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - E_trans_actin_1
     #--------------------------------------------------------------------------------------------

E_trans_actin_1_MaterialsField = iron.Field()
E_trans_actin_1_EquationsSet.MaterialsCreateStart(E_trans_actin_1_MaterialsFieldUserNumber,E_trans_actin_1_MaterialsField)

E_trans_actin_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans_actin_1 Materials Field")
E_trans_actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_actin_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
E_trans_actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_E_trans_actin)

#Diffusion Coefficient in Y
E_trans_actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_E_trans_actin)

#Storage Coefficient
E_trans_actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - E_trans_actin_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
E_trans_actin_1_SourceField = iron.Field()
E_trans_actin_1_EquationsSet.SourceCreateStart(E_trans_actin_1_SourceFieldUserNumber, E_trans_actin_1_SourceField)
E_trans_actin_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans_actin_1 Source Field")
E_trans_actin_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
E_trans_actin_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update E_trans_actin_1 Source Field
E_trans_actin_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
E_trans_actin_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update E_trans_actin_1 Dependent Field
E_trans_actin_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
E_trans_actin_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# myo_actin_1 EQUATIONS (myosin binds actin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - myo_actin_1
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
myo_actin_1_EquationsSetField         = iron.Field()
myo_actin_1_EquationsSet              = iron.EquationsSet()
myo_actin_1_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
myo_actin_1_EquationsSet.CreateStart(myo_actin_1_EquationsSetUserNumber,region,
                                   geometricField,myo_actin_1_EquationsSetSpecification,
                                   myo_actin_1_EquationsSetFieldUserNumber,myo_actin_1_EquationsSetField)
myo_actin_1_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - myo_actin_1
     #--------------------------------------------------------------------------------------------

#Start and label myo_actin_1 Field Dependent Field
myo_actin_1_Field = iron.Field()
myo_actin_1_EquationsSet.DependentCreateStart(myo_actin_1_FieldUserNumber, myo_actin_1_Field)
myo_actin_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"myo_actin_1 Field")
myo_actin_1_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"myo_actin_1 Field DELUDELN")

#Finish myo_actin_1 Field Dependent Field
myo_actin_1_EquationsSet.DependentCreateFinish()

#Initialise myo_actin_1 Field Dependent field
myo_actin_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_myo_actin_1)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - myo_actin_1
     #--------------------------------------------------------------------------------------------

myo_actin_1_MaterialsField = iron.Field()
myo_actin_1_EquationsSet.MaterialsCreateStart(myo_actin_1_MaterialsFieldUserNumber,myo_actin_1_MaterialsField)

myo_actin_1_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_actin_1 Materials Field")
myo_actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
myo_actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
myo_actin_1_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
myo_actin_1_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
myo_actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_myo_actin)

#Diffusion Coefficient in Y
myo_actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_myo_actin)

#Storage Coefficient
myo_actin_1_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)

     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - myo_actin_1
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
myo_actin_1_SourceField = iron.Field()
myo_actin_1_EquationsSet.SourceCreateStart(myo_actin_1_SourceFieldUserNumber, myo_actin_1_SourceField)
myo_actin_1_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_actin_1 Source Field")
myo_actin_1_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
myo_actin_1_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update myo_actin_1 Source Field
myo_actin_1_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
myo_actin_1_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update myo_actin_1 Dependent Field
myo_actin_1_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
myo_actin_1_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

##################################################################################################
# Membrane 2 EQUATIONS
##################################################################################################

#-------------------------------------------------------------------------------------------------
# actin_2 EQUATIONS
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - actin_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
actin_2_EquationsSetField         = iron.Field()
actin_2_EquationsSet              = iron.EquationsSet()
actin_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
actin_2_EquationsSet.CreateStart(actin_2_EquationsSetUserNumber,region,
                                   geometricField,actin_2_EquationsSetSpecification,
                                   actin_2_EquationsSetFieldUserNumber,actin_2_EquationsSetField)
actin_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - actin_2
     #--------------------------------------------------------------------------------------------

#Start and label actin_2 Field Dependent Field
actin_2_Field = iron.Field()
actin_2_EquationsSet.DependentCreateStart(actin_2_FieldUserNumber, actin_2_Field)
actin_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"actin_2 Field")
actin_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"actin_2 Field DELUDELN")

#Finish actin_2 Field Dependent Field
actin_2_EquationsSet.DependentCreateFinish()

#Initialise actin_2 Field Dependent field
actin_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_actin_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - actin_2
     #--------------------------------------------------------------------------------------------

actin_2_MaterialsField = iron.Field()
actin_2_EquationsSet.MaterialsCreateStart(actin_2_MaterialsFieldUserNumber,actin_2_MaterialsField)

actin_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_2 Materials Field")
actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
actin_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_actin)

#Diffusion Coefficient in Y
actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_actin)

#Storage Coefficient
actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - actin_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
actin_2_SourceField = iron.Field()
actin_2_EquationsSet.SourceCreateStart(actin_2_SourceFieldUserNumber, actin_2_SourceField)
actin_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_2 Source Field")
actin_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
actin_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update actin_2 Source Field
actin_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
actin_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update actin_2 Dependent Field
actin_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
actin_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# actin_a_2 EQUATIONS (Activated actin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - actin_a_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
actin_a_2_EquationsSetField         = iron.Field()
actin_a_2_EquationsSet              = iron.EquationsSet()
actin_a_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
actin_a_2_EquationsSet.CreateStart(actin_a_2_EquationsSetUserNumber,region,
                                   geometricField,actin_a_2_EquationsSetSpecification,
                                   actin_a_2_EquationsSetFieldUserNumber,actin_a_2_EquationsSetField)
actin_a_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - actin_a_2
     #--------------------------------------------------------------------------------------------

#Start and label actin_a_2 Field Dependent Field
actin_a_2_Field = iron.Field()
actin_a_2_EquationsSet.DependentCreateStart(actin_a_2_FieldUserNumber, actin_a_2_Field)
actin_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"actin_a_2 Field")
actin_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"actin_a_2 Field DELUDELN")

#Finish actin_a_2 Field Dependent Field
actin_a_2_EquationsSet.DependentCreateFinish()

#Initialise actin_a_2 Field Dependent field
actin_a_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_actin_a_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - actin_a_2
     #--------------------------------------------------------------------------------------------

actin_a_2_MaterialsField = iron.Field()
actin_a_2_EquationsSet.MaterialsCreateStart(actin_a_2_MaterialsFieldUserNumber,actin_a_2_MaterialsField)

actin_a_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_a_2 Materials Field")
actin_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
actin_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
actin_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
actin_a_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
actin_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_actin_a)

#Diffusion Coefficient in Y
actin_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_actin_a)

#Storage Coefficient
actin_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - actin_a_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
actin_a_2_SourceField = iron.Field()
actin_a_2_EquationsSet.SourceCreateStart(actin_a_2_SourceFieldUserNumber, actin_a_2_SourceField)
actin_a_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"actin_a_2 Source Field")
actin_a_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
actin_a_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update actin_a_2 Source Field
actin_a_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
actin_a_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update actin_a_2 Dependent Field
actin_a_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
actin_a_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# Rac_2 EQUATIONS
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - Rac_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
Rac_2_EquationsSetField         = iron.Field()
Rac_2_EquationsSet              = iron.EquationsSet()
Rac_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
Rac_2_EquationsSet.CreateStart(Rac_2_EquationsSetUserNumber,region,
                                   geometricField,Rac_2_EquationsSetSpecification,
                                   Rac_2_EquationsSetFieldUserNumber,Rac_2_EquationsSetField)
Rac_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - Rac_2
     #--------------------------------------------------------------------------------------------

#Start and label Rac_2 Field Dependent Field
Rac_2_Field = iron.Field()
Rac_2_EquationsSet.DependentCreateStart(Rac_2_FieldUserNumber, Rac_2_Field)
Rac_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_2 Field")
Rac_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"Rac_2 Field DELUDELN")

#Finish Rac_2 Field Dependent Field
Rac_2_EquationsSet.DependentCreateFinish()

#Initialise Rac_2 Field Dependent field
Rac_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_Rac_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - Rac_2
     #--------------------------------------------------------------------------------------------

Rac_2_MaterialsField = iron.Field()
Rac_2_EquationsSet.MaterialsCreateStart(Rac_2_MaterialsFieldUserNumber,Rac_2_MaterialsField)

Rac_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_2 Materials Field")
Rac_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
Rac_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
Rac_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
Rac_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
Rac_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_Rac)

#Diffusion Coefficient in Y
Rac_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_Rac)

#Storage Coefficient
Rac_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)

     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - Rac_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
Rac_2_SourceField = iron.Field()
Rac_2_EquationsSet.SourceCreateStart(Rac_2_SourceFieldUserNumber, Rac_2_SourceField)
Rac_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_2 Source Field")
Rac_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
Rac_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update Rac_2 Source Field
Rac_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
Rac_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update Rac_2 Dependent Field
Rac_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
Rac_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# Rac_a_2 EQUATIONS (Activated Rac)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - Rac_a_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
Rac_a_2_EquationsSetField         = iron.Field()
Rac_a_2_EquationsSet              = iron.EquationsSet()
Rac_a_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
Rac_a_2_EquationsSet.CreateStart(Rac_a_2_EquationsSetUserNumber,region,
                                   geometricField,Rac_a_2_EquationsSetSpecification,
                                   Rac_a_2_EquationsSetFieldUserNumber,Rac_a_2_EquationsSetField)
Rac_a_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - Rac_a_2
     #--------------------------------------------------------------------------------------------

#Start and label Rac_a_2 Field Dependent Field
Rac_a_2_Field = iron.Field()
Rac_a_2_EquationsSet.DependentCreateStart(Rac_a_2_FieldUserNumber, Rac_a_2_Field)
Rac_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_a_2 Field")
Rac_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"Rac_a_2 Field DELUDELN")

#Finish Rac_a_2 Field Dependent Field
Rac_a_2_EquationsSet.DependentCreateFinish()

#Initialise Rac_a_2 Field Dependent field
Rac_a_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_Rac_a_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - Rac_a_2
     #--------------------------------------------------------------------------------------------

Rac_a_2_MaterialsField = iron.Field()
Rac_a_2_EquationsSet.MaterialsCreateStart(Rac_a_2_MaterialsFieldUserNumber,Rac_a_2_MaterialsField)

Rac_a_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_a_2 Materials Field")
Rac_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
Rac_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
Rac_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
Rac_a_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
Rac_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_Rac_a)

#Diffusion Coefficient in Y
Rac_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_Rac_a)

#Storage Coefficient
Rac_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - Rac_a_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
Rac_a_2_SourceField = iron.Field()
Rac_a_2_EquationsSet.SourceCreateStart(Rac_a_2_SourceFieldUserNumber, Rac_a_2_SourceField)
Rac_a_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"Rac_a_2 Source Field")
Rac_a_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
Rac_a_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update Rac_a_2 Source Field
Rac_a_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
Rac_a_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update Rac_a_2 Dependent Field
Rac_a_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
Rac_a_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# RhoA_2 EQUATIONS (RhoA_2)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - RhoA_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
RhoA_2_EquationsSetField         = iron.Field()
RhoA_2_EquationsSet              = iron.EquationsSet()
RhoA_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
RhoA_2_EquationsSet.CreateStart(RhoA_2_EquationsSetUserNumber,region,
                                   geometricField,RhoA_2_EquationsSetSpecification,
                                   RhoA_2_EquationsSetFieldUserNumber,RhoA_2_EquationsSetField)
RhoA_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - RhoA_2
     #--------------------------------------------------------------------------------------------

#Start and label RhoA_2 Field Dependent Field
RhoA_2_Field = iron.Field()
RhoA_2_EquationsSet.DependentCreateStart(RhoA_2_FieldUserNumber, RhoA_2_Field)
RhoA_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_2 Field")
RhoA_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"RhoA_2 Field DELUDELN")

#Finish RhoA_2 Field Dependent Field
RhoA_2_EquationsSet.DependentCreateFinish()

#Initialise RhoA_2 Field Dependent field
RhoA_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_RhoA_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - RhoA_2
     #--------------------------------------------------------------------------------------------

RhoA_2_MaterialsField = iron.Field()
RhoA_2_EquationsSet.MaterialsCreateStart(RhoA_2_MaterialsFieldUserNumber,RhoA_2_MaterialsField)

RhoA_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_2 Materials Field")
RhoA_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
RhoA_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_RhoA)

#Diffusion Coefficient in Y
RhoA_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_RhoA)

#Storage Coefficient
RhoA_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - RhoA_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
RhoA_2_SourceField = iron.Field()
RhoA_2_EquationsSet.SourceCreateStart(RhoA_2_SourceFieldUserNumber, RhoA_2_SourceField)
RhoA_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_2 Source Field")
RhoA_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
RhoA_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update RhoA_2 Source Field
RhoA_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
RhoA_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update RhoA_2 Dependent Field
RhoA_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
RhoA_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)


#-------------------------------------------------------------------------------------------------
# RhoA_a_2 EQUATIONS (Activated RhoA)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - RhoA_a_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
RhoA_a_2_EquationsSetField         = iron.Field()
RhoA_a_2_EquationsSet              = iron.EquationsSet()
RhoA_a_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
RhoA_a_2_EquationsSet.CreateStart(RhoA_a_2_EquationsSetUserNumber,region,
                                   geometricField,RhoA_a_2_EquationsSetSpecification,
                                   RhoA_a_2_EquationsSetFieldUserNumber,RhoA_a_2_EquationsSetField)
RhoA_a_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - RhoA_a_2
     #--------------------------------------------------------------------------------------------

#Start and label RhoA_a_2 Field Dependent Field
RhoA_a_2_Field = iron.Field()
RhoA_a_2_EquationsSet.DependentCreateStart(RhoA_a_2_FieldUserNumber, RhoA_a_2_Field)
RhoA_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_a_2 Field")
RhoA_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"RhoA_a_2 Field DELUDELN")

#Finish RhoA_a_2 Field Dependent Field
RhoA_a_2_EquationsSet.DependentCreateFinish()

#Initialise RhoA_a_2 Field Dependent field
RhoA_a_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_RhoA_a_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - RhoA_a_2
     #--------------------------------------------------------------------------------------------

RhoA_a_2_MaterialsField = iron.Field()
RhoA_a_2_EquationsSet.MaterialsCreateStart(RhoA_a_2_MaterialsFieldUserNumber,RhoA_a_2_MaterialsField)

RhoA_a_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_a_2 Materials Field")
RhoA_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
RhoA_a_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
RhoA_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_RhoA_a)

#Diffusion Coefficient in Y
RhoA_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_RhoA_a)

#Storage Coefficient
RhoA_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - RhoA_a_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
RhoA_a_2_SourceField = iron.Field()
RhoA_a_2_EquationsSet.SourceCreateStart(RhoA_a_2_SourceFieldUserNumber, RhoA_a_2_SourceField)
RhoA_a_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"RhoA_a_2 Source Field")
RhoA_a_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
RhoA_a_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update RhoA_a_2 Source Field
RhoA_a_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
RhoA_a_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update RhoA_a_2 Dependent Field
RhoA_a_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
RhoA_a_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# myo_2 EQUATIONS (Myosin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - myo_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
myo_2_EquationsSetField         = iron.Field()
myo_2_EquationsSet              = iron.EquationsSet()
myo_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
myo_2_EquationsSet.CreateStart(myo_2_EquationsSetUserNumber,region,
                                   geometricField,myo_2_EquationsSetSpecification,
                                   myo_2_EquationsSetFieldUserNumber,myo_2_EquationsSetField)
myo_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - myo_2
     #--------------------------------------------------------------------------------------------

#Start and label myo_2 Field Dependent Field
myo_2_Field = iron.Field()
myo_2_EquationsSet.DependentCreateStart(myo_2_FieldUserNumber, myo_2_Field)
myo_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"myo_2 Field")
myo_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"myo_2 Field DELUDELN")

#Finish myo_2 Field Dependent Field
myo_2_EquationsSet.DependentCreateFinish()

#Initialise myo_2 Field Dependent field
myo_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_myo_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - myo_2
     #--------------------------------------------------------------------------------------------

myo_2_MaterialsField = iron.Field()
myo_2_EquationsSet.MaterialsCreateStart(myo_2_MaterialsFieldUserNumber,myo_2_MaterialsField)

myo_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_2 Materials Field")
myo_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
myo_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
myo_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
myo_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
myo_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_myo)

#Diffusion Coefficient in Y
myo_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_myo)

#Storage Coefficient
myo_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - myo_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
myo_2_SourceField = iron.Field()
myo_2_EquationsSet.SourceCreateStart(myo_2_SourceFieldUserNumber, myo_2_SourceField)
myo_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_2 Source Field")
myo_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
myo_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update myo_2 Source Field
myo_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
myo_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update myo_2 Dependent Field
myo_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
myo_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# myo_a_2 EQUATIONS (Activated Myosin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - myo_a_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
myo_a_2_EquationsSetField         = iron.Field()
myo_a_2_EquationsSet              = iron.EquationsSet()
myo_a_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
myo_a_2_EquationsSet.CreateStart(myo_a_2_EquationsSetUserNumber,region,
                                   geometricField,myo_a_2_EquationsSetSpecification,
                                   myo_a_2_EquationsSetFieldUserNumber,myo_a_2_EquationsSetField)
myo_a_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - myo_a_2
     #--------------------------------------------------------------------------------------------

#Start and label myo_a_2 Field Dependent Field
myo_a_2_Field = iron.Field()
myo_a_2_EquationsSet.DependentCreateStart(myo_a_2_FieldUserNumber, myo_a_2_Field)
myo_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"myo_a_2 Field")
myo_a_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"myo_a_2 Field DELUDELN")

#Finish myo_a_2 Field Dependent Field
myo_a_2_EquationsSet.DependentCreateFinish()

#Initialise myo_a_2 Field Dependent field
myo_a_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_myo_a_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - myo_a_2
     #--------------------------------------------------------------------------------------------

myo_a_2_MaterialsField = iron.Field()
myo_a_2_EquationsSet.MaterialsCreateStart(myo_a_2_MaterialsFieldUserNumber,myo_a_2_MaterialsField)

myo_a_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_a_2 Materials Field")
myo_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
myo_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
myo_a_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
myo_a_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
myo_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_myo_a)

#Diffusion Coefficient in Y
myo_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_myo_a)

#Storage Coefficient
myo_a_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - myo_a_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
myo_a_2_SourceField = iron.Field()
myo_a_2_EquationsSet.SourceCreateStart(myo_a_2_SourceFieldUserNumber, myo_a_2_SourceField)
myo_a_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_a_2 Source Field")
myo_a_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
myo_a_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update myo_a_2 Source Field
myo_a_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
myo_a_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update myo_a_2 Dependent Field
myo_a_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
myo_a_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
#  E_trans_actin_2 EQUATIONS (Bond between E-cad and activated actin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - E_trans_actin_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
E_trans_actin_2_EquationsSetField         = iron.Field()
E_trans_actin_2_EquationsSet              = iron.EquationsSet()
E_trans_actin_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
E_trans_actin_2_EquationsSet.CreateStart(E_trans_actin_2_EquationsSetUserNumber,region,
                                   geometricField,E_trans_actin_2_EquationsSetSpecification,
                                   E_trans_actin_2_EquationsSetFieldUserNumber,E_trans_actin_2_EquationsSetField)
E_trans_actin_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - E_trans_actin_2
     #--------------------------------------------------------------------------------------------

#Start and label E_trans_actin_2 Field Dependent Field
E_trans_actin_2_Field = iron.Field()
E_trans_actin_2_EquationsSet.DependentCreateStart(E_trans_actin_2_FieldUserNumber, E_trans_actin_2_Field)
E_trans_actin_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans_actin_2 Field")
E_trans_actin_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"E_trans_actin_2 Field DELUDELN")

#Finish E_trans_actin_2 Field Dependent Field
E_trans_actin_2_EquationsSet.DependentCreateFinish()

#Initialise E_trans_actin_2 Field Dependent field
E_trans_actin_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_E_trans_actin_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - E_trans_actin_2
     #--------------------------------------------------------------------------------------------

E_trans_actin_2_MaterialsField = iron.Field()
E_trans_actin_2_EquationsSet.MaterialsCreateStart(E_trans_actin_2_MaterialsFieldUserNumber,E_trans_actin_2_MaterialsField)

E_trans_actin_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans_actin_2 Materials Field")
E_trans_actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
E_trans_actin_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
E_trans_actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_E_trans_actin)

#Diffusion Coefficient in Y
E_trans_actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_E_trans_actin)

#Storage Coefficient
E_trans_actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)


     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - E_trans_actin_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
E_trans_actin_2_SourceField = iron.Field()
E_trans_actin_2_EquationsSet.SourceCreateStart(E_trans_actin_2_SourceFieldUserNumber, E_trans_actin_2_SourceField)
E_trans_actin_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"E_trans_actin_2 Source Field")
E_trans_actin_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
E_trans_actin_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update E_trans_actin_2 Source Field
E_trans_actin_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
E_trans_actin_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update E_trans_actin_2 Dependent Field
E_trans_actin_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
E_trans_actin_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# myo_actin_2 EQUATIONS (myosin binds actin)
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # EQUATION SETS - myo_actin_2
     #--------------------------------------------------------------------------------------------

# Create standard Laplace equations set
myo_actin_2_EquationsSetField         = iron.Field()
myo_actin_2_EquationsSet              = iron.EquationsSet()
myo_actin_2_EquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.REACTION_DIFFUSION_EQUATION,
                                       iron.EquationsSetSubtypes.CELLML_REAC_SPLIT_REAC_DIFF]
myo_actin_2_EquationsSet.CreateStart(myo_actin_2_EquationsSetUserNumber,region,
                                   geometricField,myo_actin_2_EquationsSetSpecification,
                                   myo_actin_2_EquationsSetFieldUserNumber,myo_actin_2_EquationsSetField)
myo_actin_2_EquationsSet.CreateFinish()



     #--------------------------------------------------------------------------------------------
     # DEPENDENT FIELD - myo_actin_2
     #--------------------------------------------------------------------------------------------

#Start and label myo_actin_2 Field Dependent Field
myo_actin_2_Field = iron.Field()
myo_actin_2_EquationsSet.DependentCreateStart(myo_actin_2_FieldUserNumber, myo_actin_2_Field)
myo_actin_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"myo_actin_2 Field")
myo_actin_2_Field.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"myo_actin_2 Field DELUDELN")

#Finish myo_actin_2 Field Dependent Field
myo_actin_2_EquationsSet.DependentCreateFinish()

#Initialise myo_actin_2 Field Dependent field
myo_actin_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,init_myo_actin_2)


     #--------------------------------------------------------------------------------------------
     # MATERIAL FIELD - myo_actin_2
     #--------------------------------------------------------------------------------------------

myo_actin_2_MaterialsField = iron.Field()
myo_actin_2_EquationsSet.MaterialsCreateStart(myo_actin_2_MaterialsFieldUserNumber,myo_actin_2_MaterialsField)

myo_actin_2_MaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_actin_2 Materials Field")
myo_actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.NODE_BASED)
myo_actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.NODE_BASED)
myo_actin_2_MaterialsField.ComponentInterpolationSet(iron.FieldVariableTypes.U,3,iron.FieldInterpolationTypes.NODE_BASED)
myo_actin_2_EquationsSet.MaterialsCreateFinish()


#Diffusion Coefficient in X
myo_actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     1,Dx_myo_actin)

#Diffusion Coefficient in Y
myo_actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     2,Dy_myo_actin)

#Storage Coefficient
myo_actin_2_MaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                     iron.FieldParameterSetTypes.VALUES,
                                                     3,store_coeff)

     #--------------------------------------------------------------------------------------------
     # SOURCE FIELD - myo_actin_2
     #--------------------------------------------------------------------------------------------

#Setting up the source field for reaction diffusion equation set.
#For split problem subtype, source field is not used
myo_actin_2_SourceField = iron.Field()
myo_actin_2_EquationsSet.SourceCreateStart(myo_actin_2_SourceFieldUserNumber, myo_actin_2_SourceField)
myo_actin_2_SourceField.VariableLabelSet(iron.FieldVariableTypes.U,"myo_actin_2 Source Field")
myo_actin_2_EquationsSet.SourceCreateFinish()

#Initialise to 0.0
myo_actin_2_SourceField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                  iron.FieldParameterSetTypes.VALUES,
                                                  1,0.0)


#Update myo_actin_2 Source Field
myo_actin_2_SourceField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                              iron.FieldParameterSetTypes.VALUES)
myo_actin_2_SourceField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES)

#Update myo_actin_2 Dependent Field
myo_actin_2_Field.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES)
myo_actin_2_Field.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                         iron.FieldParameterSetTypes.VALUES)

#-------------------------------------------------------------------------------------------------
# CONSTANT FIELDS
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # k0_1 FIELD
     #--------------------------------------------------------------------------------------------

k0_1_Field = iron.Field()
k0_1_Field.CreateStart(k0_1_FieldUserNumber,region)

k0_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k0_1_Field.MeshDecompositionSet(decomposition)
k0_1_Field.GeometricFieldSet(geometricField)
k0_1_Field.NumberOfVariablesSet(1)
k0_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k0_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k0_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k0_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k0_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k0_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k0_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k0_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k0_1_Field.CreateFinish()

k0_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k0_1)

     #--------------------------------------------------------------------------------------------
     # k1_1 FIELD
     #--------------------------------------------------------------------------------------------

k1_1_Field = iron.Field()
k1_1_Field.CreateStart(k1_1_FieldUserNumber,region)

k1_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k1_1_Field.MeshDecompositionSet(decomposition)
k1_1_Field.GeometricFieldSet(geometricField)
k1_1_Field.NumberOfVariablesSet(1)
k1_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k1_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k1_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k1_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k1_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k1_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k1_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k1_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k1_1_Field.CreateFinish()

k1_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k1_1)

     #--------------------------------------------------------------------------------------------
     # k2_1 FIELD
     #--------------------------------------------------------------------------------------------

k2_1_Field = iron.Field()
k2_1_Field.CreateStart(k2_1_FieldUserNumber,region)

k2_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k2_1_Field.MeshDecompositionSet(decomposition)
k2_1_Field.GeometricFieldSet(geometricField)
k2_1_Field.NumberOfVariablesSet(1)
k2_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k2_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k2_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k2_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k2_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k2_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k2_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k2_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k2_1_Field.CreateFinish()

k2_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k2_1)

     #--------------------------------------------------------------------------------------------
     # k3_1 FIELD
     #--------------------------------------------------------------------------------------------

k3_1_Field = iron.Field()
k3_1_Field.CreateStart(k3_1_FieldUserNumber,region)

k3_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k3_1_Field.MeshDecompositionSet(decomposition)
k3_1_Field.GeometricFieldSet(geometricField)
k3_1_Field.NumberOfVariablesSet(1)
k3_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k3_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k3_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k3_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k3_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k3_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k3_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k3_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k3_1_Field.CreateFinish()

k3_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k3_1)

     #--------------------------------------------------------------------------------------------
     # k3r_1 FIELD
     #--------------------------------------------------------------------------------------------

k3r_1_Field = iron.Field()
k3r_1_Field.CreateStart(k3r_1_FieldUserNumber,region)

k3r_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k3r_1_Field.MeshDecompositionSet(decomposition)
k3r_1_Field.GeometricFieldSet(geometricField)
k3r_1_Field.NumberOfVariablesSet(1)
k3r_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k3r_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k3r_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k3r_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k3r_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k3r_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k3r_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k3r_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k3r_1_Field.CreateFinish()

k3r_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k3r_1)



     #--------------------------------------------------------------------------------------------
     # k4_1 FIELD
     #--------------------------------------------------------------------------------------------

k4_1_Field = iron.Field()
k4_1_Field.CreateStart(k4_1_FieldUserNumber,region)

k4_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k4_1_Field.MeshDecompositionSet(decomposition)
k4_1_Field.GeometricFieldSet(geometricField)
k4_1_Field.NumberOfVariablesSet(1)
k4_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k4_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k4_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k4_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k4_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k4_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k4_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k4_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k4_1_Field.CreateFinish()

k4_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k4_1)

     #--------------------------------------------------------------------------------------------
     # k5_1 FIELD
     #--------------------------------------------------------------------------------------------

k5_1_Field = iron.Field()
k5_1_Field.CreateStart(k5_1_FieldUserNumber,region)

k5_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k5_1_Field.MeshDecompositionSet(decomposition)
k5_1_Field.GeometricFieldSet(geometricField)
k5_1_Field.NumberOfVariablesSet(1)
k5_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k5_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k5_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k5_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k5_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k5_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k5_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k5_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k5_1_Field.CreateFinish()

k5_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k5_1)
 
     #--------------------------------------------------------------------------------------------
     # k5r_1 FIELD
     #--------------------------------------------------------------------------------------------

k5r_1_Field = iron.Field()
k5r_1_Field.CreateStart(k5r_1_FieldUserNumber,region)

k5r_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k5r_1_Field.MeshDecompositionSet(decomposition)
k5r_1_Field.GeometricFieldSet(geometricField)
k5r_1_Field.NumberOfVariablesSet(1)
k5r_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k5r_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k5r_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k5r_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k5r_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k5r_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k5r_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k5r_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k5r_1_Field.CreateFinish()

k5r_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k5r_1)


     #--------------------------------------------------------------------------------------------
     # k6_1 FIELD
     #--------------------------------------------------------------------------------------------

k6_1_Field = iron.Field()
k6_1_Field.CreateStart(k6_1_FieldUserNumber,region)

k6_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k6_1_Field.MeshDecompositionSet(decomposition)
k6_1_Field.GeometricFieldSet(geometricField)
k6_1_Field.NumberOfVariablesSet(1)
k6_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k6_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k6_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k6_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k6_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k6_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k6_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k6_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k6_1_Field.CreateFinish()

k6_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k6_1)

     #--------------------------------------------------------------------------------------------
     # k7_1 FIELD
     #--------------------------------------------------------------------------------------------

k7_1_Field = iron.Field()
k7_1_Field.CreateStart(k7_1_FieldUserNumber,region)

k7_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k7_1_Field.MeshDecompositionSet(decomposition)
k7_1_Field.GeometricFieldSet(geometricField)
k7_1_Field.NumberOfVariablesSet(1)
k7_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k7_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k7_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k7_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k7_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k7_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k7_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k7_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k7_1_Field.CreateFinish()

k7_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k7_1)

     #--------------------------------------------------------------------------------------------
     # k7r_1 FIELD
     #--------------------------------------------------------------------------------------------

k7r_1_Field = iron.Field()
k7r_1_Field.CreateStart(k7r_1_FieldUserNumber,region)

k7r_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k7r_1_Field.MeshDecompositionSet(decomposition)
k7r_1_Field.GeometricFieldSet(geometricField)
k7r_1_Field.NumberOfVariablesSet(1)
k7r_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k7r_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k7r_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k7r_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k7r_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k7r_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k7r_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k7r_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k7r_1_Field.CreateFinish()

k7r_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k7r_1)

     #--------------------------------------------------------------------------------------------
     # k8_1 FIELD
     #--------------------------------------------------------------------------------------------

k8_1_Field = iron.Field()
k8_1_Field.CreateStart(k8_1_FieldUserNumber,region)

k8_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k8_1_Field.MeshDecompositionSet(decomposition)
k8_1_Field.GeometricFieldSet(geometricField)
k8_1_Field.NumberOfVariablesSet(1)
k8_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k8_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k8_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k8_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k8_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k8_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k8_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k8_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k8_1_Field.CreateFinish()

k8_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k8_1)

     #--------------------------------------------------------------------------------------------
     # k9_1 FIELD
     #--------------------------------------------------------------------------------------------

k9_1_Field = iron.Field()
k9_1_Field.CreateStart(k9_1_FieldUserNumber,region)

k9_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k9_1_Field.MeshDecompositionSet(decomposition)
k9_1_Field.GeometricFieldSet(geometricField)
k9_1_Field.NumberOfVariablesSet(1)
k9_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k9_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k9_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k9_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k9_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k9_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k9_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k9_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k9_1_Field.CreateFinish()

k9_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k9_1)

     #--------------------------------------------------------------------------------------------
     # k9r_1 FIELD
     #--------------------------------------------------------------------------------------------

k9r_1_Field = iron.Field()
k9r_1_Field.CreateStart(k9r_1_FieldUserNumber,region)

k9r_1_Field.TypeSet(iron.FieldTypes.GENERAL)
k9r_1_Field.MeshDecompositionSet(decomposition)
k9r_1_Field.GeometricFieldSet(geometricField)
k9r_1_Field.NumberOfVariablesSet(1)
k9r_1_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k9r_1_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k9r_1_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k9r_1_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k9r_1_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k9r_1 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k9r_1_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k9r_1_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k9r_1_Field.CreateFinish()

k9r_1_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k9r_1)

##################################################################################################
# Membrane 2
##################################################################################################

     #--------------------------------------------------------------------------------------------
     # k0_2 FIELD
     #--------------------------------------------------------------------------------------------

k0_2_Field = iron.Field()
k0_2_Field.CreateStart(k0_2_FieldUserNumber,region)

k0_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k0_2_Field.MeshDecompositionSet(decomposition)
k0_2_Field.GeometricFieldSet(geometricField)
k0_2_Field.NumberOfVariablesSet(1)
k0_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k0_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k0_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k0_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k0_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k0_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k0_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k0_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k0_2_Field.CreateFinish()

k0_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k0_2)

     #--------------------------------------------------------------------------------------------
     # k1_2 FIELD
     #--------------------------------------------------------------------------------------------

k1_2_Field = iron.Field()
k1_2_Field.CreateStart(k1_2_FieldUserNumber,region)

k1_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k1_2_Field.MeshDecompositionSet(decomposition)
k1_2_Field.GeometricFieldSet(geometricField)
k1_2_Field.NumberOfVariablesSet(1)
k1_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k1_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k1_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k1_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k1_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k1_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k1_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k1_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k1_2_Field.CreateFinish()

k1_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k1_2)

     #--------------------------------------------------------------------------------------------
     # k2_2 FIELD
     #--------------------------------------------------------------------------------------------

k2_2_Field = iron.Field()
k2_2_Field.CreateStart(k2_2_FieldUserNumber,region)

k2_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k2_2_Field.MeshDecompositionSet(decomposition)
k2_2_Field.GeometricFieldSet(geometricField)
k2_2_Field.NumberOfVariablesSet(1)
k2_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k2_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k2_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k2_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k2_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k2_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k2_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k2_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k2_2_Field.CreateFinish()

k2_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k2_2)

     #--------------------------------------------------------------------------------------------
     # k3_2 FIELD
     #--------------------------------------------------------------------------------------------

k3_2_Field = iron.Field()
k3_2_Field.CreateStart(k3_2_FieldUserNumber,region)

k3_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k3_2_Field.MeshDecompositionSet(decomposition)
k3_2_Field.GeometricFieldSet(geometricField)
k3_2_Field.NumberOfVariablesSet(1)
k3_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k3_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k3_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k3_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k3_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k3_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k3_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k3_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k3_2_Field.CreateFinish()

k3_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k3_2)

     #--------------------------------------------------------------------------------------------
     # k3r_2 FIELD
     #--------------------------------------------------------------------------------------------

k3r_2_Field = iron.Field()
k3r_2_Field.CreateStart(k3r_2_FieldUserNumber,region)

k3r_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k3r_2_Field.MeshDecompositionSet(decomposition)
k3r_2_Field.GeometricFieldSet(geometricField)
k3r_2_Field.NumberOfVariablesSet(1)
k3r_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k3r_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k3r_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k3r_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k3r_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k3r_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k3r_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k3r_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k3r_2_Field.CreateFinish()

k3r_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k3r_2)



     #--------------------------------------------------------------------------------------------
     # k4_2 FIELD
     #--------------------------------------------------------------------------------------------

k4_2_Field = iron.Field()
k4_2_Field.CreateStart(k4_2_FieldUserNumber,region)

k4_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k4_2_Field.MeshDecompositionSet(decomposition)
k4_2_Field.GeometricFieldSet(geometricField)
k4_2_Field.NumberOfVariablesSet(1)
k4_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k4_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k4_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k4_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k4_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k4_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k4_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k4_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k4_2_Field.CreateFinish()

k4_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k4_2)

     #--------------------------------------------------------------------------------------------
     # k5_2 FIELD
     #--------------------------------------------------------------------------------------------

k5_2_Field = iron.Field()
k5_2_Field.CreateStart(k5_2_FieldUserNumber,region)

k5_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k5_2_Field.MeshDecompositionSet(decomposition)
k5_2_Field.GeometricFieldSet(geometricField)
k5_2_Field.NumberOfVariablesSet(1)
k5_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k5_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k5_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k5_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k5_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k5_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k5_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k5_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k5_2_Field.CreateFinish()

k5_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k5_2)
 
     #--------------------------------------------------------------------------------------------
     # k5r_2 FIELD
     #--------------------------------------------------------------------------------------------

k5r_2_Field = iron.Field()
k5r_2_Field.CreateStart(k5r_2_FieldUserNumber,region)

k5r_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k5r_2_Field.MeshDecompositionSet(decomposition)
k5r_2_Field.GeometricFieldSet(geometricField)
k5r_2_Field.NumberOfVariablesSet(1)
k5r_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k5r_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k5r_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k5r_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k5r_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k5r_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k5r_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k5r_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k5r_2_Field.CreateFinish()

k5r_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k5r_2)


     #--------------------------------------------------------------------------------------------
     # k6_2 FIELD
     #--------------------------------------------------------------------------------------------

k6_2_Field = iron.Field()
k6_2_Field.CreateStart(k6_2_FieldUserNumber,region)

k6_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k6_2_Field.MeshDecompositionSet(decomposition)
k6_2_Field.GeometricFieldSet(geometricField)
k6_2_Field.NumberOfVariablesSet(1)
k6_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k6_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k6_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k6_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k6_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k6_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k6_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k6_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k6_2_Field.CreateFinish()

k6_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k6_2)

     #--------------------------------------------------------------------------------------------
     # k7_2 FIELD
     #--------------------------------------------------------------------------------------------

k7_2_Field = iron.Field()
k7_2_Field.CreateStart(k7_2_FieldUserNumber,region)

k7_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k7_2_Field.MeshDecompositionSet(decomposition)
k7_2_Field.GeometricFieldSet(geometricField)
k7_2_Field.NumberOfVariablesSet(1)
k7_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k7_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k7_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k7_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k7_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k7_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k7_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k7_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k7_2_Field.CreateFinish()

k7_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k7_2)

     #--------------------------------------------------------------------------------------------
     # k7r_2 FIELD
     #--------------------------------------------------------------------------------------------

k7r_2_Field = iron.Field()
k7r_2_Field.CreateStart(k7r_2_FieldUserNumber,region)

k7r_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k7r_2_Field.MeshDecompositionSet(decomposition)
k7r_2_Field.GeometricFieldSet(geometricField)
k7r_2_Field.NumberOfVariablesSet(1)
k7r_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k7r_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k7r_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k7r_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k7r_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k7r_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k7r_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k7r_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k7r_2_Field.CreateFinish()

k7r_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k7r_2)

     #--------------------------------------------------------------------------------------------
     # k8_2 FIELD
     #--------------------------------------------------------------------------------------------

k8_2_Field = iron.Field()
k8_2_Field.CreateStart(k8_2_FieldUserNumber,region)

k8_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k8_2_Field.MeshDecompositionSet(decomposition)
k8_2_Field.GeometricFieldSet(geometricField)
k8_2_Field.NumberOfVariablesSet(1)
k8_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k8_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k8_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k8_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k8_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k8_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k8_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k8_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k8_2_Field.CreateFinish()

k8_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k8_2)

     #--------------------------------------------------------------------------------------------
     # k9_2 FIELD
     #--------------------------------------------------------------------------------------------

k9_2_Field = iron.Field()
k9_2_Field.CreateStart(k9_2_FieldUserNumber,region)

k9_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k9_2_Field.MeshDecompositionSet(decomposition)
k9_2_Field.GeometricFieldSet(geometricField)
k9_2_Field.NumberOfVariablesSet(1)
k9_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k9_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k9_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k9_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k9_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k9_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k9_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k9_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k9_2_Field.CreateFinish()

k9_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k9_2)

     #--------------------------------------------------------------------------------------------
     # k9r_2 FIELD
     #--------------------------------------------------------------------------------------------

k9r_2_Field = iron.Field()
k9r_2_Field.CreateStart(k9r_2_FieldUserNumber,region)

k9r_2_Field.TypeSet(iron.FieldTypes.GENERAL)
k9r_2_Field.MeshDecompositionSet(decomposition)
k9r_2_Field.GeometricFieldSet(geometricField)
k9r_2_Field.NumberOfVariablesSet(1)
k9r_2_Field.VariableTypesSet([iron.FieldVariableTypes.U])
k9r_2_Field.DataTypeSet(iron.FieldVariableTypes.U,iron.FieldDataTypes.DP)
k9r_2_Field.DimensionSet(iron.FieldVariableTypes.U,iron.FieldDimensionTypes.SCALAR)
k9r_2_Field.NumberOfComponentsSet(iron.FieldVariableTypes.U,1)
k9r_2_Field.VariableLabelSet(iron.FieldVariableTypes.U,"k9r_2 Field")

geometricMeshComponent = geometricField.ComponentMeshComponentGet(iron.FieldVariableTypes.U,1)
k9r_2_Field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,geometricMeshComponent)
k9r_2_Field.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,
                                      iron.FieldInterpolationTypes.NODE_BASED)

k9r_2_Field.CreateFinish()

k9r_2_Field.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,
                                        1,k9r_2)

#-------------------------------------------------------------------------------------------------
# CELLML
#-------------------------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------------------
     # CELLML FIELD
     #--------------------------------------------------------------------------------------------

#Initialise cellml
cellml = iron.CellML()
cellml.CreateStart(cellmlUserNumber,region)


#Importing the cellml model
constantModelIndex = cellml.ModelImport(Actin_Ecad)

# Kinetic Constants
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k0_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k1_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k2_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k3_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k3r_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k4_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k5_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k5r_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k6_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k7_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k7r_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k8_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k9_1")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k9r_1")

cellml.VariableSetAsKnown(constantModelIndex,"membrane/k0_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k1_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k2_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k3_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k3r_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k4_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k5_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k5r_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k6_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k7_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k7r_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k8_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k9_2")
cellml.VariableSetAsKnown(constantModelIndex,"membrane/k9r_2")

#Finish cellml
cellml.CreateFinish()

     #--------------------------------------------------------------------------------------------
     # CELLML OPENCMISS FIELD MAPS
     #--------------------------------------------------------------------------------------------

#Initialise Field Maps
cellml.FieldMapsCreateStart()

#------------------------#
# Kinetic Constants #
#------------------------#

#k0_1
cellml.CreateFieldToCellMLMap(k0_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k0_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k0_1",
                              iron.FieldParameterSetTypes.VALUES,k0_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k1_1
cellml.CreateFieldToCellMLMap(k1_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k1_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k1_1",
                              iron.FieldParameterSetTypes.VALUES,k1_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k2_1
cellml.CreateFieldToCellMLMap(k2_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k2_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k2_1",
                              iron.FieldParameterSetTypes.VALUES,k2_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#k3_1
cellml.CreateFieldToCellMLMap(k3_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k3_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k3_1",
                              iron.FieldParameterSetTypes.VALUES,k3_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k3r_1
cellml.CreateFieldToCellMLMap(k3r_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k3r_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k3r_1",
                              iron.FieldParameterSetTypes.VALUES,k3r_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k4_1
cellml.CreateFieldToCellMLMap(k4_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k4_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k4_1",
                              iron.FieldParameterSetTypes.VALUES,k4_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#k5_1
cellml.CreateFieldToCellMLMap(k5_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k5_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k5_1",
                              iron.FieldParameterSetTypes.VALUES,k5_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#k5r_1
cellml.CreateFieldToCellMLMap(k5r_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k5r_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k5r_1",
                              iron.FieldParameterSetTypes.VALUES,k5r_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k6_1
cellml.CreateFieldToCellMLMap(k6_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k6_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k6_1",
                              iron.FieldParameterSetTypes.VALUES,k6_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k7_1
cellml.CreateFieldToCellMLMap(k7_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k7_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k7_1",
                              iron.FieldParameterSetTypes.VALUES,k7_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k7r_1
cellml.CreateFieldToCellMLMap(k7r_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k7r_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k7r_1",
                              iron.FieldParameterSetTypes.VALUES,k7r_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k8_1
cellml.CreateFieldToCellMLMap(k8_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k8_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k8_1",
                              iron.FieldParameterSetTypes.VALUES,k8_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
    
#k9_1
cellml.CreateFieldToCellMLMap(k9_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k9_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k9_1",
                              iron.FieldParameterSetTypes.VALUES,k9_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k9r_1
cellml.CreateFieldToCellMLMap(k9r_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k9r_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k9r_1",
                              iron.FieldParameterSetTypes.VALUES,k9r_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k0_2
cellml.CreateFieldToCellMLMap(k0_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k0_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k0_2",
                              iron.FieldParameterSetTypes.VALUES,k0_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k1_2
cellml.CreateFieldToCellMLMap(k1_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k1_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k1_2",
                              iron.FieldParameterSetTypes.VALUES,k1_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k2_2
cellml.CreateFieldToCellMLMap(k2_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k2_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k2_2",
                              iron.FieldParameterSetTypes.VALUES,k2_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#k3_2
cellml.CreateFieldToCellMLMap(k3_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k3_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k3_2",
                              iron.FieldParameterSetTypes.VALUES,k3_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k3r_2
cellml.CreateFieldToCellMLMap(k3r_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k3r_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k3r_2",
                              iron.FieldParameterSetTypes.VALUES,k3r_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k4_2
cellml.CreateFieldToCellMLMap(k4_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k4_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k4_2",
                              iron.FieldParameterSetTypes.VALUES,k4_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#k5_2
cellml.CreateFieldToCellMLMap(k5_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k5_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k5_2",
                              iron.FieldParameterSetTypes.VALUES,k5_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#k5r_2
cellml.CreateFieldToCellMLMap(k5r_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k5r_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k5r_2",
                              iron.FieldParameterSetTypes.VALUES,k5r_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k6_2
cellml.CreateFieldToCellMLMap(k6_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k6_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k6_2",
                              iron.FieldParameterSetTypes.VALUES,k6_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k7_2
cellml.CreateFieldToCellMLMap(k7_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k7_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k7_2",
                              iron.FieldParameterSetTypes.VALUES,k7_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k7r_2
cellml.CreateFieldToCellMLMap(k7r_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k7r_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k7r_2",
                              iron.FieldParameterSetTypes.VALUES,k7r_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
#k8_2
cellml.CreateFieldToCellMLMap(k8_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k8_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k8_2",
                              iron.FieldParameterSetTypes.VALUES,k8_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
    
#k9_2
cellml.CreateFieldToCellMLMap(k9_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k9_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k9_2",
                              iron.FieldParameterSetTypes.VALUES,k9_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#k9r_2
cellml.CreateFieldToCellMLMap(k9r_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/k9r_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/k9r_2",
                              iron.FieldParameterSetTypes.VALUES,k9r_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)


#---------------------#
# Dependent Variables #
#---------------------#

#E_trans
cellml.CreateFieldToCellMLMap(E_trans_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/E_trans",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/E_trans",
                              iron.FieldParameterSetTypes.VALUES,E_trans_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

#actin_1
cellml.CreateFieldToCellMLMap(actin_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/actin_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/actin_1",
                              iron.FieldParameterSetTypes.VALUES,actin_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#actin_a_1
cellml.CreateFieldToCellMLMap(actin_a_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/actin_a_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/actin_a_1",
                              iron.FieldParameterSetTypes.VALUES,actin_a_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             
    
#Rac_1
cellml.CreateFieldToCellMLMap(Rac_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/Rac_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/Rac_1",
                              iron.FieldParameterSetTypes.VALUES,Rac_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#Rac_a_1
cellml.CreateFieldToCellMLMap(Rac_a_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/Rac_a_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/Rac_a_1",
                              iron.FieldParameterSetTypes.VALUES,Rac_a_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#RhoA_1
cellml.CreateFieldToCellMLMap(RhoA_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/RhoA_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/RhoA_1",
                              iron.FieldParameterSetTypes.VALUES,RhoA_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#RhoA_a_1
cellml.CreateFieldToCellMLMap(RhoA_a_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/RhoA_a_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/RhoA_a_1",
                              iron.FieldParameterSetTypes.VALUES,RhoA_a_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#myo_1
cellml.CreateFieldToCellMLMap(myo_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/myo_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/myo_1",
                              iron.FieldParameterSetTypes.VALUES,myo_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#myo_a_1
cellml.CreateFieldToCellMLMap(myo_a_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/myo_a_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/myo_a_1",
                              iron.FieldParameterSetTypes.VALUES,myo_a_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#E_trans_actin_1
cellml.CreateFieldToCellMLMap(E_trans_actin_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/E_trans_actin_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/E_trans_actin_1",
                              iron.FieldParameterSetTypes.VALUES,E_trans_actin_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#myo_actin_1
cellml.CreateFieldToCellMLMap(myo_actin_1_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/myo_actin_1",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/myo_actin_1",
                              iron.FieldParameterSetTypes.VALUES,myo_actin_1_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#actin_2
cellml.CreateFieldToCellMLMap(actin_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/actin_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/actin_2",
                              iron.FieldParameterSetTypes.VALUES,actin_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
                              
#actin_a_2
cellml.CreateFieldToCellMLMap(actin_a_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/actin_a_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/actin_a_2",
                              iron.FieldParameterSetTypes.VALUES,actin_a_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             
    
#Rac_2
cellml.CreateFieldToCellMLMap(Rac_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/Rac_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/Rac_2",
                              iron.FieldParameterSetTypes.VALUES,Rac_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#Rac_a_2
cellml.CreateFieldToCellMLMap(Rac_a_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/Rac_a_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/Rac_a_2",
                              iron.FieldParameterSetTypes.VALUES,Rac_a_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#RhoA_2
cellml.CreateFieldToCellMLMap(RhoA_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/RhoA_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/RhoA_2",
                              iron.FieldParameterSetTypes.VALUES,RhoA_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#RhoA_a_2
cellml.CreateFieldToCellMLMap(RhoA_a_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/RhoA_a_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/RhoA_a_2",
                              iron.FieldParameterSetTypes.VALUES,RhoA_a_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#myo_2
cellml.CreateFieldToCellMLMap(myo_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/myo_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/myo_2",
                              iron.FieldParameterSetTypes.VALUES,myo_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#myo_a_2
cellml.CreateFieldToCellMLMap(myo_a_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/myo_a_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/myo_a_2",
                              iron.FieldParameterSetTypes.VALUES,myo_a_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#E_trans_actin_2
cellml.CreateFieldToCellMLMap(E_trans_actin_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/E_trans_actin_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/E_trans_actin_2",
                              iron.FieldParameterSetTypes.VALUES,E_trans_actin_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)                             

#myo_actin_2
cellml.CreateFieldToCellMLMap(myo_actin_2_Field,iron.FieldVariableTypes.U,1,
                              iron.FieldParameterSetTypes.VALUES,constantModelIndex,
                              "membrane/myo_actin_2",iron.FieldParameterSetTypes.VALUES)
cellml.CreateCellMLToFieldMap(constantModelIndex,"membrane/myo_actin_2",
                              iron.FieldParameterSetTypes.VALUES,myo_actin_2_Field,
                              iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)   
#Finish Field Maps
cellml.FieldMapsCreateFinish()

     #--------------------------------------------------------------------------------------------
     # MODEL FIELD - CELLML
     #--------------------------------------------------------------------------------------------

#Initialise the Model Field
cellmlModelsField = iron.Field()
cellml.ModelsFieldCreateStart(cellmlModelsFieldUserNumber,cellmlModelsField)
cellml.ModelsFieldCreateFinish()


     #--------------------------------------------------------------------------------------------
     # STATE FIELD - CELLML
     #--------------------------------------------------------------------------------------------

#Initialise the State Field
cellmlStateField = iron.Field()
cellml.StateFieldCreateStart(cellmlStateFieldUserNumber,cellmlStateField)
cellml.StateFieldCreateFinish()

     #-------------------------------------------------------------------------------------------------
     # PARAMETERS FIELD - CELLML
     #-------------------------------------------------------------------------------------------------

#Initialise the Parameters Field
cellmlParametersField = iron.Field()
cellml.ParametersFieldCreateStart(cellmlParametersFieldUserNumber,cellmlParametersField)
cellml.ParametersFieldCreateFinish()

#-------------------------------------------------------------------------------------------------
# REACTION DIFFUSION EQUATIONS
#-------------------------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------
    # E_trans
    #-------------------------------------------------------------------------------------------------

#E_trans Equations Set
E_trans_Equations = iron.Equations()
E_trans_EquationsSet.EquationsCreateStart(E_trans_Equations)
E_trans_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
E_trans_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
E_trans_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # actin_1
    #-------------------------------------------------------------------------------------------------

#actin_1 Equations Set
actin_1_Equations = iron.Equations()
actin_1_EquationsSet.EquationsCreateStart(actin_1_Equations)
actin_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
actin_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
actin_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # actin_a_1
    #-------------------------------------------------------------------------------------------------

#actin_a_1 Equations Set
actin_a_1_Equations = iron.Equations()
actin_a_1_EquationsSet.EquationsCreateStart(actin_a_1_Equations)
actin_a_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
actin_a_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
actin_a_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # Rac_1
    #-------------------------------------------------------------------------------------------------

#Rac_1 Equations Set
Rac_1_Equations = iron.Equations()
Rac_1_EquationsSet.EquationsCreateStart(Rac_1_Equations)
Rac_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Rac_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
Rac_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # Rac_a_1
    #-------------------------------------------------------------------------------------------------

#Rac_a_1 Equations Set
Rac_a_1_Equations = iron.Equations()
Rac_a_1_EquationsSet.EquationsCreateStart(Rac_a_1_Equations)
Rac_a_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Rac_a_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
Rac_a_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # RhoA_1
    #-------------------------------------------------------------------------------------------------

#RhoA_1 Equations Set
RhoA_1_Equations = iron.Equations()
RhoA_1_EquationsSet.EquationsCreateStart(RhoA_1_Equations)
RhoA_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
RhoA_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
RhoA_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # RhoA_a_1
    #-------------------------------------------------------------------------------------------------

#RhoA_a_1 Equations Set
RhoA_a_1_Equations = iron.Equations()
RhoA_a_1_EquationsSet.EquationsCreateStart(RhoA_a_1_Equations)
RhoA_a_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
RhoA_a_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
RhoA_a_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # myo_1
    #-------------------------------------------------------------------------------------------------

#myo_1 Equations Set
myo_1_Equations = iron.Equations()
myo_1_EquationsSet.EquationsCreateStart(myo_1_Equations)
myo_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
myo_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
myo_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # myo_a_1
    #-------------------------------------------------------------------------------------------------

#myo_a_1 Equations Set
myo_a_1_Equations = iron.Equations()
myo_a_1_EquationsSet.EquationsCreateStart(myo_a_1_Equations)
myo_a_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
myo_a_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
myo_a_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # E_trans_actin_1
    #-------------------------------------------------------------------------------------------------

#E_trans_actin_1 Equations Set
E_trans_actin_1_Equations = iron.Equations()
E_trans_actin_1_EquationsSet.EquationsCreateStart(E_trans_actin_1_Equations)
E_trans_actin_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
E_trans_actin_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
E_trans_actin_1_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # myo_actin_1
    #-------------------------------------------------------------------------------------------------

#myo_actin_1 Equations Set
myo_actin_1_Equations = iron.Equations()
myo_actin_1_EquationsSet.EquationsCreateStart(myo_actin_1_Equations)
myo_actin_1_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
myo_actin_1_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
myo_actin_1_EquationsSet.EquationsCreateFinish()


    #-------------------------------------------------------------------------------------------------
    # actin_2
    #-------------------------------------------------------------------------------------------------

#actin_2 Equations Set
actin_2_Equations = iron.Equations()
actin_2_EquationsSet.EquationsCreateStart(actin_2_Equations)
actin_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
actin_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
actin_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # actin_a_2
    #-------------------------------------------------------------------------------------------------

#actin_a_2 Equations Set
actin_a_2_Equations = iron.Equations()
actin_a_2_EquationsSet.EquationsCreateStart(actin_a_2_Equations)
actin_a_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
actin_a_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
actin_a_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # Rac_2
    #-------------------------------------------------------------------------------------------------

#Rac_2 Equations Set
Rac_2_Equations = iron.Equations()
Rac_2_EquationsSet.EquationsCreateStart(Rac_2_Equations)
Rac_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Rac_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
Rac_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # Rac_a_2
    #-------------------------------------------------------------------------------------------------

#Rac_a_2 Equations Set
Rac_a_2_Equations = iron.Equations()
Rac_a_2_EquationsSet.EquationsCreateStart(Rac_a_2_Equations)
Rac_a_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
Rac_a_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
Rac_a_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # RhoA_2
    #-------------------------------------------------------------------------------------------------

#RhoA_2 Equations Set
RhoA_2_Equations = iron.Equations()
RhoA_2_EquationsSet.EquationsCreateStart(RhoA_2_Equations)
RhoA_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
RhoA_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
RhoA_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # RhoA_a_2
    #-------------------------------------------------------------------------------------------------

#RhoA_a_2 Equations Set
RhoA_a_2_Equations = iron.Equations()
RhoA_a_2_EquationsSet.EquationsCreateStart(RhoA_a_2_Equations)
RhoA_a_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
RhoA_a_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
RhoA_a_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # myo_2
    #-------------------------------------------------------------------------------------------------

#myo_2 Equations Set
myo_2_Equations = iron.Equations()
myo_2_EquationsSet.EquationsCreateStart(myo_2_Equations)
myo_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
myo_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
myo_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # myo_a_2
    #-------------------------------------------------------------------------------------------------

#myo_a_2 Equations Set
myo_a_2_Equations = iron.Equations()
myo_a_2_EquationsSet.EquationsCreateStart(myo_a_2_Equations)
myo_a_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
myo_a_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
myo_a_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # E_trans_actin_2
    #-------------------------------------------------------------------------------------------------

#E_trans_actin_2 Equations Set
E_trans_actin_2_Equations = iron.Equations()
E_trans_actin_2_EquationsSet.EquationsCreateStart(E_trans_actin_2_Equations)
E_trans_actin_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
E_trans_actin_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
E_trans_actin_2_EquationsSet.EquationsCreateFinish()

    #-------------------------------------------------------------------------------------------------
    # myo_actin_2
    #-------------------------------------------------------------------------------------------------

#myo_actin_2 Equations Set
myo_actin_2_Equations = iron.Equations()
myo_actin_2_EquationsSet.EquationsCreateStart(myo_actin_2_Equations)
myo_actin_2_Equations.sparsityType = iron.EquationsSparsityTypes.SPARSE

#Equations Outputs
myo_actin_2_Equations.outputType = iron.EquationsOutputTypes.NONE

#Finish
myo_actin_2_EquationsSet.EquationsCreateFinish()


#_________________________________________________________________________________________________
# REACTION DIFFUSION PROBLEM AND SOLVER
#_________________________________________________________________________________________________

#-------------------------------------------------------------------------------------------------
# PROBLEM - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------
    
# Create Problem
rd_problem = iron.Problem()
rd_problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
                           iron.ProblemTypes.REACTION_DIFFUSION_EQUATION,
                           iron.ProblemSubtypes.CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT]
rd_problem.CreateStart(rd_problemUserNumber, rd_problemSpecification)
rd_problem.CreateFinish()


# Create control loops
rd_problem.ControlLoopCreateStart()
rd_controlLoop = iron.ControlLoop()

rd_problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],rd_controlLoop)
# start time, stop time, time increment
rd_controlLoop.TimesSet(startT,startT+deltat*freq_update,deltat*freq_update)


#Control Loop Outputs
#rd_controlLoop.TimeOutputSet(outputfreq)
#rd_controlLoop.LoadOutputSet(1)
#rd_controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
#rd_controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)

rd_problem.ControlLoopCreateFinish()


#-------------------------------------------------------------------------------------------------
# SOLVER - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

#
#    1st Solver --> DAE
#         |
#         v
#    2nd Solver --> Dynamic 
#         |
#         v      
#    3rd Solver --> DAE
#

#Create rd_problem solver for Strang splitting
rd_problem.SolversCreateStart()


#Create first solver --> DAE Solver
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,rd_solver)
#rd_solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
rd_solver.DAETimeStepSet(ODE_Timestep)
rd_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Create second solver --> Dynamic solver for parabolic equation
rd_solver = iron.Solver()
rd_linearsolver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,rd_solver)

#Set theta - backward vs forward time step parameter
rd_solver.DynamicThetaSet([1.0])

#Set output type
rd_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Obtain dynamic linear solver from the solver
rd_solver.DynamicLinearSolverGet(rd_linearsolver)

#Set Library
#rd_solver.LibraryTypeSet(iron.SolverLibraries.LAPACK)
#rd_solver.LibraryTypeSet(iron.SolverLibraries.CMISS)
#rd_solver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
#rd_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)

#rd_solver.LinearDirectTypeSet(iron.DirectLinearSolverTypes.LU)

rd_linearsolver.LinearIterativeMaximumIterationsSet(100000)
#rd_linearsolver.linearIterativeAbsoluteTolerance = 1.0E-12
#rd_linearsolver.linearIterativeRelativeTolerance = 1.0E-12


#Create third solver --> Another DAE Solver
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,rd_solver)
#rd_solver.DAESolverTypeSet(iron.DAESolverTypes.EULER)
rd_solver.DAETimeStepSet(ODE_Timestep)
rd_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)

#Finish the rd_problem
rd_problem.SolversCreateFinish()

#-------------------------------------------------------------------------------------------------
# SOLVER CELLML EQUATIONS - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

#Start rd_solver Cellml Equations
rd_problem.CellMLEquationsCreateStart()

#Create first solver
#cellml equations
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,rd_solver)
cellmlEquations = iron.CellMLEquations()
rd_solver.CellMLEquationsGet(cellmlEquations)
#Add into Cellml Environment
cellmlIndex = cellmlEquations.CellMLAdd(cellml)

#Create third solver
#cellml equations
rd_solver = iron.Solver()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,rd_solver)
cellmlEquations = iron.CellMLEquations()
rd_solver.CellMLEquationsGet(cellmlEquations)
#Add into Cellml Environment
cellmlIndex = cellmlEquations.CellMLAdd(cellml)

#Finish the solver Cellml Equations
rd_problem.CellMLEquationsCreateFinish()

#-------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

#Start solver equations
rd_problem.SolverEquationsCreateStart()

#Create second solver
#Solver equations
rd_solver = iron.Solver()
rd_solverEquations = iron.SolverEquations()
rd_problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,rd_solver)
rd_solver.SolverEquationsGet(rd_solverEquations)


rd_solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE

equationsSetIndex = rd_solverEquations.EquationsSetAdd(E_trans_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(actin_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(actin_a_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(Rac_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(Rac_a_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(RhoA_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(RhoA_a_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(myo_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(myo_a_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(E_trans_actin_1_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(myo_actin_1_EquationsSet)

equationsSetIndex = rd_solverEquations.EquationsSetAdd(actin_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(actin_a_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(Rac_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(Rac_a_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(RhoA_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(RhoA_a_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(myo_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(myo_a_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(E_trans_actin_2_EquationsSet)
equationsSetIndex = rd_solverEquations.EquationsSetAdd(myo_actin_2_EquationsSet)
rd_problem.SolverEquationsCreateFinish()


#-------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS - REACTION DIFFUSION
#-------------------------------------------------------------------------------------------------

# Set up Boundary Conditions & Perturbation Initial Conditions
rd_boundaryConditions = iron.BoundaryConditions()
rd_solverEquations.BoundaryConditionsCreateStart(rd_boundaryConditions)
rd_solverEquations.BoundaryConditionsCreateFinish()


#-------------------------------------------------------------------------------------------------
# End of OPENCMISS SETTINGS 
#------------------------------------------------------------------------------------------------- 

#-------------------------------------------------------------------------------------------------
# START OF LATTICE-BASED MODEL & COUPLING BETWEEN TWO MODELS 
#------------------------------------------------------------------------------------------------- 


#----------------------------------------------------------------------------------
# Correspondence between lattices and the elements
#----------------------------------------------------------------------------------
Area_all = np.loadtxt('AreaAll.txt')

matchedlist = []
with open('MatchLatticeElement_7.0.txt','r') as matchedfile:
    for line in matchedfile:
        matchedlist.append(list(map(int,line.rstrip().split('   '))))

# All the lattice indexes on the mesh
contact_lattice = []
for items in matchedlist:
    for item in items:
        contact_lattice.append(item)
        
contact_lattice_index = np.unique(np.array(contact_lattice))
contact_lattice_TF = np.zeros((x_len**2,1),dtype=int)
contact_lattice_TF[contact_lattice_index] = 1

#----------------------------------------------------------------------------------
# NODE INDEX ON THE EDGE & CENTER
#----------------------------------------------------------------------------------

#Find and calculate nodes' index & the number of node on the edge and in the central area 
node_edge, node_edge_left, node_edge_right, node_central =  (np.array([],dtype=int),
                                                            np.array([],dtype=int), 
                                                            np.array([],dtype=int), 
                                                            np.array([],dtype=int))
node_edge_yesorno = np.zeros(num_nodes_incenter,dtype=int)

for i in range(num_nodes_incenter):
    node_rad = np.sqrt(node_incenter[i+1,1]**2+node_incenter[i+1,2]**2)       # Distances between nodes and the center
    if node_rad > cortex_rad:
        node_edge = np.append(node_edge,i+1)
        if node_incenter[i+1,1] >=0:
            node_edge_right = np.append(node_edge_right,len(node_edge)-1)
        else: 
            node_edge_left = np.append(node_edge_left,len(node_edge)-1)
    else:
        node_central = np.append(node_central,i+1)

node_edge_yesorno[node_edge-1] = 1 # set all index of nodes on the edge to be one
num_node_edge = node_edge.size
num_node_central = node_central.size




#----------------------------------------------------------------------------------
# NEIGHBOR NODES
#----------------------------------------------------------------------------------
# Find neighbor nodes of a node,
# based on the angle formed by the F-actin (3um, 3.0/5.0)

theta_threshold = float(Fila_len)/real_radius   # F-actin effective region expressed as an angle formed by a F-actin
theta_reverse = np.arccos(node_incenter[node_edge,1]/np.sqrt(np.sum(node_incenter[node_edge,1:3]**2,1))) # theta (radial) location of each node on the edge 
theta_reverse[node_incenter[node_edge,2]<0] = 2*np.pi-theta_reverse[node_incenter[node_edge,2]<0]   # for y < 0, change their theta from 0<theta<pi to pi<theta<2pi

theta_nearby_node = []

node_fixed_angle = np.zeros(num_node_edge,dtype=float) 

node_fixed_angle[np.logical_or((theta_reverse<np.pi/2),(theta_reverse>np.pi/2*3))] = 1.0

Top90angleNode = [] # node index
Bot90angleNode = [] # node index 
Top90angle = 0  # average contact angle of node between 90 and 100 degrees, unit: rad
Bot90angle = 0  # average contact angle of node between 260 and 270 degrees, unit: rad

# Find nearby nodes of the node within F-actin effective region
for i in range(num_node_edge):
    # Find nearby node of the current node
    theta_diff = np.abs(theta_reverse-theta_reverse[i])
    nearby_node = node_edge[theta_diff<(theta_threshold/2.0)]
    nearby_node_index = np.array(np.where(theta_diff<(theta_threshold/2.0)))
    theta_nearby_node.append(node_edge[nearby_node_index])
    
    # Find the node of Top90angle and the node of Bot90angle
    if theta_reverse[i]>(np.pi/2.0) and theta_reverse[i]<(np.pi/2.0+np.deg2rad(10)):
        Top90angleNode.append(i)
    elif theta_reverse[i]>(np.pi/2.0+np.pi-np.deg2rad(10)) and theta_reverse[i]<(np.pi/2.0+np.pi):
        Bot90angleNode.append(i)

lattice_nearby_node = [] # corresponding lattice points on nearby nodes 

for i in range(num_node_edge):        
    
    # count all lattice points on the nearby elements
    all_lattice = np.array([],dtype = int)
    for j in range(len(theta_nearby_node[i].ravel())):
        all_lattice = np.append(all_lattice,matchedlist[theta_nearby_node[i].ravel()[j]-1])
    lattice_nearby_node.append(all_lattice)    

#----------------------------------------------------------------------------------
# NORMAL VECTOR OF ELEMENTS (VECTORS PARALLEL TO THE RIM)
#----------------------------------------------------------------------------------
# Get the vectors parallel to the rim of the contact at each element on the edge. 

# Clockwise vectors
vectors_clockwise = np.array([[np.cos(rad+np.pi/2.0),np.sin(rad+np.pi/2.0)] for rad in theta_reverse])
# Counterclockwise vectors
vectors_reverselyclockwise = np.array([[np.cos(rad-np.pi/2.0),np.sin(rad-np.pi/2.0)] for rad in theta_reverse])

        
#----------------------------------------------------------------------------------
# E-cad species initialization Section
#----------------------------------------------------------------------------------


total_side_len = side_len*x_len         # Side length of the square
radius_nm = total_side_len/2.0             # radius of the circle in nm
num_sites = x_len * y_len                   # The number of lattices on the surface

gamma_t = math.exp(g_trans)
ass_p_t = gamma_t / (1 + gamma_t)   # Association probability trans
dis_p_t = 1.0 / (1 + gamma_t)         # Dissociation probability trans

gamma_t_ringright = math.exp(g_trans_ringright)
ass_p_t_ringright = gamma_t_ringright / (1 + gamma_t_ringright)   # Association probability trans
dis_p_t_ringright = 1.0 / (1 + gamma_t_ringright)         # Dissociation probability trans

gamma_c = math.exp(g_cis)
ass_p_c = gamma_c / (1 + gamma_c)   # Association probability cis
dis_p_c = 1.0 / (1 + gamma_c)         # Dissociation probability cis

#Initialize position and orientations of M1 and M2
# col_0: index, 
# col_1: monomer, 
# col_2: trans, 
# col_3: cis donor
# col_4: orientation(N: 5, E: 6, S: 7, W: 8),
# col_5: x, 
# col_6: y, 
# col_7: occupied,
# col_8: cis receiver
# col_9: E_trans_actin bond (bound == 1, free == 0)

M1 = [[0 for x in range(num_sites)] for y in range(10)]  # Initialize M1 matrix
M2 = [[0 for x in range(num_sites)] for y in range(10)]  # Initialize M2 matrix

# assign index and x,y positions for lattices
for i in range(x_len*y_len):
    M1[0][i] = i
    M1[5][i] = i % x_len
    M1[6][i] = int(math.floor(i/x_len))

M2[0][:], M2[5][:], M2[6][:] = M1[0][:], M1[5][:], M1[6][:]

# Find all lattice points within the circle
par_x, par_y = [0 for x in range(num_sites)],[0 for x in range(num_sites)]
if x_len%2 != 0:
    par_x = (np.array(M1[5][:])-x_len/2)*side_len
    par_y = (np.array(M1[6][:])-y_len/2)*side_len
else:
    par_x = (np.array(M1[5][:])-(x_len/2-0.5))*side_len
    par_y = (np.array(M1[6][:])-(y_len/2-0.5))*side_len
dist_center = np.sqrt(par_x**2+par_y**2)

contact_lattice_index = np.unique(np.array(contact_lattice))
contact_lattice_TF = np.zeros((x_len**2,1),dtype=int)
contact_lattice_TF[contact_lattice_index] = 1

# edit properties of occupied sites
num_par_M1_st = int(Ecad_conc_M1*x_len**2)
num_par_M2_st = int(Ecad_conc_M2*x_len**2)

sites_index = list(range(x_len*y_len))
occupied_site_M1 = [select(sites_index) for x in range(num_par_M1_st)]  # randomly select num_par occupied sites,M1
sites_index = list(range(x_len*y_len))
occupied_site_M2 = [select(sites_index) for x in range(num_par_M2_st)]  # randomly select num_par occupied sites,M2

for i in range(num_par_M1_st):
    M1[1][occupied_site_M1[i]] = 1
    M1[4][occupied_site_M1[i]] = random.randint(5, 8)
    M1[7][occupied_site_M1[i]] = 1
    
for i in range(num_par_M2_st):    
    M2[1][occupied_site_M2[i]] = 1
    M2[4][occupied_site_M2[i]] = random.randint(5, 8)
    M2[7][occupied_site_M2[i]] = 1

M1 = np.transpose(np.array(M1))
M2 = np.transpose(np.array(M2))

onedge = M1[np.logical_and((dist_center<radius_nm),(dist_center>cortex_rad*1000)),0]
incenter = M1[(dist_center<cortex_rad*1000),0]
area_onedge = np.pi*((radius_um)**2-(cortex_rad)**2)
area_incenter = np.pi*cortex_rad**2

col_index = np.array([1, 2, 3, 4, 7, 8, 9])


###############################################################################
#  Whole contact 
###############################################################################

incont = np.array(np.where((np.sqrt(par_x**2+par_y**2)<(radius_um)*1000)))        # Indexes of lattice points within the circle
incont_list = list(map(int,incont.ravel().tolist()))                  # turn it into a list
incont_TF = np.zeros(x_len*y_len,dtype=int)
incont_TF[incont] = 1                                         # Label lattices inside the circle by 1

###############################################################################
#  Central contact 
###############################################################################

incent = np.array(np.where((np.sqrt(par_x**2+par_y**2)<(radius_um-cort_thick)*1000)))        # Indexes of lattice points within the circle
incent_list = list(map(int,incent.ravel().tolist()))                  # turn it into a list
incent_TF = np.zeros(x_len*y_len,dtype=int)
incent_TF[incent] = 1                                         # Label lattices inside the circle by 1

###############################################################################
#  Left Half Ring 
###############################################################################

inring_left = np.array(np.where((np.sqrt(par_x**2+par_y**2)<radius_um*1000)&(np.sqrt(par_x**2+par_y**2)>(radius_um-cort_thick)*1000)&(par_x<0)))        # Indexes of lattice points within the circle
inring_left_list = list(map(int,inring_left.ravel().tolist()))                  # turn it into a list
inring_left_TF = np.zeros(x_len*y_len,dtype=int)
inring_left_TF[inring_left] = 1                                         # Label lattices inside the circle by 1

###############################################################################
#  Right Half Ring 
###############################################################################

inring_right = np.array(np.where((np.sqrt(par_x**2+par_y**2)<radius_um*1000)&(np.sqrt(par_x**2+par_y**2)>(radius_um-cort_thick)*1000)&(par_x>=0)))        # Indexes of lattice points within the circle
inring_right_list = list(map(int,inring_right.ravel().tolist()))                  # turn it into a list
inring_right_TF = np.zeros(x_len*y_len,dtype=int)
inring_right_TF[inring_right] = 1                                         # Label lattices inside the circle by 1


#-------------------------------------------------------------------------------------------------
# TIME LOOP FOR THE E-CAD DYNAMICS
#-------------------------------------------------------------------------------------------------


file_M1 = open("M1_data.txt", "w")
file_M2 = open("M2_data.txt", "w")
file_con = open("con_all.txt", "w")
file_num = open("num_trans_cis_immobile.txt", "w")
num_cad_all = np.zeros((int(t_total_step/record_freq+1),19))
edgeangle = []
currentT = startT

t_100 = 0

for t in range(t_total_step+1):
    # top-bottom or left-right
    incre = [1,x_len]
    # top or bottom, left or right
    oneminusone = [1,-1]
    
    # Filter occupied site to simplify the calculation   
    M1_occupied = M1[np.array(np.where(M1[:,7]==1)).ravel(),:]
    M2_occupied = M2[np.array(np.where(M2[:,7]==1)).ravel(),:]
    
    #######################################################################
    # Formation of dissociation of cis and trans
    # (1: monomer, 2: trans, 3: Cis, 4: direction)
    #######################################################################
    
    num_par_M1 = np.count_nonzero(M1[:,7])
    num_par_M2 = np.count_nonzero(M2[:,7])

    # Formation of trans dimer
    # monomer to trans (1 to 2)
    
    meet_mono_trans = np.array(np.where((M1[:,7]==1)&(M2[:,7]==1)&(M1[:,2]!=1)&((M1[:,4]-M2[:,4]==3)|(M1[:,4]-M2[:,4]==-1)))).ravel()

    new_trans = meet_mono_trans[np.array([random.random() for x in range(meet_mono_trans.size)]) < ass_p_t]
    new_trans = new_trans[(contact_lattice_TF[new_trans]!=0).ravel()]   # Exclude all trans dimers outside the contact area

    # Formation of Cis

    # initialize two arrays
    cis_pos_M1_ass = np.zeros(num_par_M1,dtype=int)
    cis_pos_M2_ass = np.zeros(num_par_M2,dtype=int)
    # loop through all occupied lattice, check posible positions for cis association 
    # top-bottom or left-right
    for i in range(num_par_M1):
        cis_pos_M1_ass[i] = M1_occupied[i,0]+incre[M1_occupied[i,4]%2]*oneminusone[M1_occupied[i,4]//7]
    for i in range(num_par_M2):
        cis_pos_M2_ass[i] = M2_occupied[i,0]+incre[M2_occupied[i,4]%2]*oneminusone[M2_occupied[i,4]//7]
    
    # periodic boundary condition
    cis_pos_M1_ass[(M1_occupied[:,0]%x_len==(x_len-1)) & (cis_pos_M1_ass%x_len==0)] -= x_len
    cis_pos_M1_ass[(M1_occupied[:,0]%x_len==0) & (cis_pos_M1_ass%x_len==(x_len-1))] += x_len
    cis_pos_M1_ass[cis_pos_M1_ass<0] += x_len*y_len 
    cis_pos_M1_ass[cis_pos_M1_ass>x_len*y_len-1] -= x_len*y_len 
    
    cis_pos_M2_ass[(M2_occupied[:,0]%x_len==(x_len-1)) & (cis_pos_M2_ass%x_len==0)] -= x_len
    cis_pos_M2_ass[(M2_occupied[:,0]%x_len==0) & (cis_pos_M2_ass%x_len==(x_len-1))] += x_len
    cis_pos_M2_ass[cis_pos_M2_ass<0] += x_len*y_len 
    cis_pos_M2_ass[cis_pos_M2_ass>x_len*y_len-1] -= x_len*y_len
    
    # if the direction of E-cad on the chosen position doesn't equal to the direction 
    # of E-cad of choosing position, unselect the chose position.        
    for i in range(num_par_M1):
        if M1[cis_pos_M1_ass[i],4] != M1_occupied[i,4]:
            cis_pos_M1_ass[i] = 0
    for i in range(num_par_M2):
        if M2[cis_pos_M2_ass[i],4] != M2_occupied[i,4]:
            cis_pos_M2_ass[i] = 0
    
    # find index of lattices that already form the cis
    exist_cis_M1 = np.array(np.where(M1_occupied[:,3]==1)).ravel()
    exist_cis_M2 = np.array(np.where(M2_occupied[:,3]==1)).ravel()
    
    # eliminate exist cis bonds from cis chosen list
    cis_pos_M1_tobeass = np.copy(cis_pos_M1_ass)
    cis_pos_M1_tobeass[exist_cis_M1] = 0
    
    cis_pos_M2_tobeass = np.copy(cis_pos_M2_ass)
    cis_pos_M2_tobeass[exist_cis_M2] = 0
    
    # set all points which doesn't pass MC to be zeros
    # Cis Association index
    cis_pos_M1_tobeass[np.array([random.random() for x in range(num_par_M1)]) > ass_p_c] = 0  
    cis_pos_M2_tobeass[np.array([random.random() for x in range(num_par_M2)]) > ass_p_c] = 0
    
    # Number of cis bonds association
    num_cis_M1_ass = np.count_nonzero(cis_pos_M1_tobeass)
    num_cis_M2_ass = np.count_nonzero(cis_pos_M2_tobeass)
    
    # Cis Dissociation index
    cis_dis_index_M1 = exist_cis_M1[np.array([random.random() for x in range(exist_cis_M1.size)]) < dis_p_c]
    cis_dis_index_M2 = exist_cis_M2[np.array([random.random() for x in range(exist_cis_M2.size)]) < dis_p_c]
    
    # Number of cis bonds dissociation
    num_cis_M1_dis = cis_dis_index_M1.size
    num_cis_M2_dis = cis_dis_index_M2.size
    
    # Trans Dissociation

    trans_index = np.array(np.where(M1[:,2]==1)).ravel()
    trans_dis = trans_index[np.array([random.random() for x in range(trans_index.size)]) < dis_p_t]     # index of dissociated trans
    

    ###########################################################################
    # update positions states 'moving direction' orientations 
    # monomer: position states 'moving direction' orientation
    # trans: position states 'moving direction' 'both orientation'
    # cis: states
    ###########################################################################

    # Cis Dissociation
    M1[M1_occupied[cis_dis_index_M1,0],3] = 0
    M1[cis_pos_M1_ass[cis_dis_index_M1],8] = 0
    
    M2[M2_occupied[cis_dis_index_M2,0],3] = 0
    M2[cis_pos_M2_ass[cis_dis_index_M2],8] = 0
    
    # Cis Association
    M1[M1_occupied[cis_pos_M1_tobeass!=0,0],3] = 1         # choosing site
    M1[cis_pos_M1_tobeass[cis_pos_M1_tobeass!=0],8] = 1    # choosing site
    
    M2[M2_occupied[cis_pos_M2_tobeass!=0,0],3] = 1         # choosing site
    M2[cis_pos_M2_tobeass[cis_pos_M2_tobeass!=0],8] = 1    # choosing site
    
    # update trans formation
    M1[new_trans,2] = 1
    M2[new_trans,2] = 1
    
    # update trans dissociation
    M1[trans_dis,2] = 0
    M2[trans_dis,2] = 0

    # update monomer
    # Dissociated trans or cis. if a site is occupied but is not the state is not trans or cis, it must be a mono
    M1[(M1[:,2]==0)&(M1[:,3]==0)&(M1[:,8]==0)&(M1[:,7]==1),1] = 1
    M2[(M2[:,2]==0)&(M2[:,3]==0)&(M2[:,8]==0)&(M2[:,7]==1),1] = 1
    # Associated Monomers. if the state of a site is trans or cis, it must not be a mono
    M1[(M1[:,3]!=0)|(M1[:,8]!=0)|(M1[:,2]!=0),1] = 0
    M2[(M2[:,3]!=0)|(M2[:,8]!=0)|(M2[:,2]!=0),1] = 0
    
    #----------------------------------------------------------------------
    # update M1 position
    #----------------------------------------------------------------------
    M1_occupied_new = np.array(np.where(M1[:,7]==1))
    mono_index_M1 = np.transpose(np.array(np.where((M1[:,1]==1)&(M1[:,9]==0)))).ravel()    # index of free monomers
    moving_dir_M1 = [random.randint(5, 8) for x in range(mono_index_M1.size)]
    
    M1_pos = np.c_[mono_index_M1-1, mono_index_M1+1, mono_index_M1-x_len, mono_index_M1+x_len]
    
    # check for periodic boundary conditions
    M1_pos[M1_pos[:,0]%x_len==(x_len-1),0] += x_len
    M1_pos[M1_pos[:,1]%x_len==0,1] -= x_len
    M1_pos[M1_pos[:,2]<0,2] += x_len*y_len
    M1_pos[M1_pos[:,3]>x_len*y_len-1,3] -= x_len*y_len
    
    #----------------------------------------------------------------------
    # update M2 position
    #----------------------------------------------------------------------
    M2_occupied_new = np.array(np.where(M2[:,7]==1))
    mono_index_M2 = np.transpose(np.array(np.where((M2[:,1]==1)&(M2[:,9]==0)))).ravel()     # index of free monomers
    moving_dir_M2 = [random.randint(5, 8) for x in range(mono_index_M2.size)]
    
    M2_pos = np.c_[mono_index_M2-1, mono_index_M2+1, mono_index_M2-x_len, mono_index_M2+x_len]
    
    # check for periodic boundary conditions
    M2_pos[M2_pos[:,0]%x_len==(x_len-1),0] += x_len
    M2_pos[M2_pos[:,1]%x_len==0,1] -= x_len
    M2_pos[M2_pos[:,2]<0,2] += x_len*y_len
    M2_pos[M2_pos[:,3]>x_len*y_len-1,3] -= x_len*y_len

    #----------------------------------------------------------------------
    # update trans position
    #----------------------------------------------------------------------
    all_occupied_new = np.unique(np.append(M1_occupied_new, M2_occupied_new))       # all occupied sites on M1 and M2
    trans_index_movable = np.transpose(np.array(np.where((M1[:,2]==1)&(M2[:,2]==1)&(M1[:,3]==0)&(M2[:,3]==0)&(M1[:,8]==0)&(M2[:,8]==0)&(M1[:,9]==0)&(M2[:,9]==0)))).ravel()     # Trans dimers that are not in cis or not in actin binding
    moving_dir_trans = [random.randint(5, 8) for x in range(trans_index_movable.size)]
    
    trans_pos = np.c_[trans_index_movable-1, trans_index_movable+1, trans_index_movable-x_len, trans_index_movable+x_len]

    # check for periodic boundary conditions
    trans_pos[trans_pos[:,0]%x_len==(x_len-1),0] += x_len
    trans_pos[trans_pos[:,1]%x_len==0,1] -= x_len
    trans_pos[trans_pos[:,2]<0,2] += x_len*y_len
    trans_pos[trans_pos[:,3]>x_len*y_len-1,3] -= x_len*y_len
    
    #----------------------------------------------------------------------
    # update all randomly
    #----------------------------------------------------------------------
    trans_pos[contact_lattice_TF.ravel()[trans_pos]==0] = 0
    num_M1 = mono_index_M1.size
    num_M2 = mono_index_M2.size
    num_trans_temp = trans_index_movable.size
    random_select = list(range(num_M1+num_M2+num_trans_temp))
    
    for n in range(num_M1+num_M2+num_trans_temp):
        m = select(random_select)
        if m<num_M1:
            i = m
            # check if lattice points around the monomer are available and if they are in the circle
            
            conditions = M1[M1_pos[i,:],7] == 0

            conditions[M1_pos[i,:]==0] = False
            
            random_move = random.randint(0, 3)
            desti = M1_pos[i, random_move]     # randint is inclusive at both ends
            
            # Reject the move if move to the occupied site or out of the circle
            if conditions[random_move] == True:
                M1[desti,col_index] = M1[mono_index_M1[i],col_index]
                M1[mono_index_M1[i],col_index] = 0
            else: 
                continue
                
        elif m>num_M1-1 and m<num_M1+num_M2:
            i=m-num_M1
        # check if lattice points around the monomer are available and if they are in the circle
            conditions = M2[M2_pos[i,:],7] == 0
            conditions[M2_pos[i,:]==0] = False
            
            random_move = random.randint(0, 3)
            desti = M2_pos[i, random_move]     # randint is inclusive at both ends
            
            # Reject the move if move to the occupied site or out of the circle
            if conditions[random_move] == True:
                M2[desti,col_index] = M2[mono_index_M2[i],col_index]
                M2[mono_index_M2[i],col_index] = 0
            else: 
                continue
             
        elif m<num_M1+num_M2+num_trans_temp and m > num_M1+num_M2-1:
            i = m-num_M1-num_M2
            
            conditions = np.logical_and(M1[trans_pos[i,:],7] == 0,M2[trans_pos[i,:],7] == 0)
            conditions[trans_pos[i,:]==0] = False
            
            random_move = random.randint(0, 3)
            desti = trans_pos[i, random_move]     # randint is inclusive at both ends
            
            # Reject the move if move to the occupied site or out of the circle
            if conditions[random_move] == True:
                M1[desti,col_index] = M1[trans_index_movable[i],col_index]
                M1[trans_index_movable[i],col_index] = 0
                M2[desti,col_index] = M2[trans_index_movable[i],col_index]
                M2[trans_index_movable[i],col_index] = 0
            else: 
                continue
    
    M1[(M1[:,1]==1)&(M1[:,9]==0),4] = [random.randint(5, 8) for x in range(mono_index_M1.size)]          # update monomer orientation
    M2[(M2[:,1]==1)&(M2[:,9]==0),4] = [random.randint(5, 8) for x in range(mono_index_M2.size)]          # update monomer orientation
    
    trans_index_movable = np.transpose(np.array(np.where((M1[:,2]==1)&(M2[:,2]==1)&(M1[:,3]==0)&(M2[:,3]==0)&(M1[:,8]==0)&(M2[:,8]==0)&(M1[:,9]==0)&(M2[:,9]==0)))).ravel()     # Trans dimers that are not in cis

    M1[trans_index_movable,4] = [random.randint(5, 8) for x in range(trans_index_movable.size)]          # update orientation
    M2[trans_index_movable,4] = M1[trans_index_movable,4]+1     # The direction of M2 equals to the direction of M1 rotate -90 degrees
    M2[M2[:,4]==9,4] = M2[M2[:,4]==9,4]-4                       # replace direction 4 by 8
    
    for i in col_index:
        M1[M1[:,7] == 0, i] = 0
        M2[M2[:,7] == 0, i] = 0
    
    
    # top-bottom or left-right
    incre = [1,x_len]
    # top or bottom, left or right
    oneminusone = [1,-1]
    
    # Filter occupied site to simplify the calculation   
    M1_occupied = M1[np.array(np.where(M1[:,7]==1)).ravel(),:]
    M2_occupied = M2[np.array(np.where(M2[:,7]==1)).ravel(),:]

    # Record M1 and M2 in txt
    if t%record_freq==0:
        file_M1.write('{} {} {} {} {} {} {} {} {} {}\n'.format(M1_occupied[:,1].size,0,0,0,0,0,0,0,0,0))
        file_M2.write('{} {} {} {} {} {} {} {} {} {}\n'.format(M2_occupied[:,1].size,0,0,0,0,0,0,0,0,0))
        for i in range(M1_occupied[:,1].size):
            string_M1 =  "    ".join(np.around(M1_occupied[i,:],decimals=2).astype(str))+"\n"
            file_M1.write(string_M1)
            
        for i in range(M2_occupied[:,1].size):
            string_M2 =  "    ".join(np.around(M2_occupied[i,:],decimals=2).astype(str))+"\n"
            file_M2.write(string_M2)
        
    # -------------------------------------------------------------------------------------------------
    # COUPLED the Reaction-diffusion model with the lattice based model
    # -------------------------------------------------------------------------------------------------

    if t%freq_update == 0:

        print("Time: %i"%t)
        currentT = currentT+deltat*freq_update
        
        SF = radius_um/real_radius
        All_F_on_Ecad = F_cort*0.2985
        
        # Number of E-cad before (prev timestep)
        num_bound_onedge_before = np.array(np.where((M1[onedge,2]==1)&(M1[onedge,9]!=0))).size
        num_bound_incenter_before = np.array(np.where((M1[incenter,2]==1)&(M1[incenter,9]!=0))).size
        
    	  # initialize variables
        broken_bond = np.array([])
        new_bond = np.array([])
        np.zeros(num_nodes_incenter,dtype=float)
        N_trans_unbound = np.zeros(num_nodes_incenter,dtype=float)
        # N_E_trans_actin_continuum = []
        N_E_trans_actin_lattice_1 = 0
        N_E_trans_actin_lattice_2 = 0
        
        num_broken_bond_1 = np.zeros(num_nodes_incenter,dtype=float)
        num_broken_bond_2 = np.zeros(num_nodes_incenter,dtype=float)
        num_new_bond = np.zeros(num_nodes_incenter,dtype=float)
        N_actin_a_all = np.zeros(num_nodes_incenter,dtype=float)
        N_E_trans_actin_1 = np.zeros(num_nodes_incenter,dtype=int)
        N_E_trans_actin_2 = np.zeros(num_nodes_incenter,dtype=int)
        
        N_trans_unbound_all = np.zeros(num_nodes_incenter,dtype=float)
        N_new_bond_all = np.zeros(num_nodes_incenter,dtype=float)
        F_alpha_all = np.array([])
        # Scale factors calculated based on the myosin concentration
        k4_SF_edge = np.array([])
        k4_SF_center = np.array([])
        k4_SF_tilt = np.array([])
        k4_SF_untilt = np.array([])

        k10r_all = []
        
        con_myo_tilt = np.array([])

#        n_prob_all = np.zeros(num_nodes_incenter,dtype=float)
        # Calculate the Concentration of E-cad trans dimers in each element
        E_trans_con_all = []
  
        # -------------------------------------------------------------------------------------------------
        # Calculate force from the central region & dissociation of E-cad-trans-actin bond on the central region
        # -------------------------------------------------------------------------------------------------
        
        num_actin_a_1 = 0.0
        num_myo_actin_1 = 0.0
        
        con_myo_actin_1 = np.zeros(num_nodes_incenter,dtype=float)    # Initialization 
        con_actin_a_1 = np.zeros(num_nodes_incenter,dtype=float)      # Initialization
        con_myo_actin_2 = np.zeros(num_nodes_incenter,dtype=float)    # Initialization 
        con_actin_a_2 = np.zeros(num_nodes_incenter,dtype=float)      # Initialization
        
        con_E_trans = np.zeros(num_nodes_incenter,dtype=float)    # Initialization 
        
        for i in range(num_nodes_incenter):
            # -------------------------------------------------------------------------------------------------
            # Retrive concentration of current element
            # -------------------------------------------------------------------------------------------------
            
            con_actin_a_1[i] = actin_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
            con_myo_actin_1[i] = myo_actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
            con_actin_a_2[i] = actin_a_2_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
            con_myo_actin_2[i] = myo_actin_2_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
            
            con_E_trans[i] = E_trans_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
            
            num_actin_a_1 = num_actin_a_1+con_actin_a_1[i]*Area_all[i]
            num_myo_actin_1 = num_myo_actin_1+con_myo_actin_1[i]*Area_all[i]
            
            # index of lattice point in the ith (current) element
            indexoflattice = np.array(matchedlist[i]) 

            # -------------------------------------------------------------------------------------------------
            # Formation of E-cad/Actin on the whole field
            # -------------------------------------------------------------------------------------------------
            form_p_1 = deltat*freq_update*k5_1*con_actin_a_1[i]   # Probability of E_cad_actin bond formation
            form_p_2 = deltat*freq_update*k5_2*con_actin_a_2[i]   # Probability of E_cad_actin bond formation
            
            # exist E-cad                    
            exist_cad_1 = indexoflattice[M1[indexoflattice,9]==0]
            exist_cad_2 = indexoflattice[M2[indexoflattice,9]==0]
            
            N_exist_cad_1 = exist_cad_1.size
            N_exist_cad_2 = exist_cad_2.size
            
            # exist E-cad-actin of the branched network on lattice field      
            exist_bond_1 = indexoflattice[M1[indexoflattice,9]==1]
            exist_bond_2 = indexoflattice[M2[indexoflattice,9]==1]
            N_exist_bond_1 = exist_bond_1.size
            N_exist_bond_2 = exist_bond_2.size
            
            # New Ecad/F-actin Bonds Formation 
            new_bonds_random_1 = np.random.rand(N_exist_cad_1)
            new_bonds_1 = exist_cad_1[new_bonds_random_1<form_p_1]
            M1[new_bonds_1,9] = 1
            
            new_bonds_random_2 = np.random.rand(N_exist_cad_2)
            new_bonds_2 = exist_cad_2[new_bonds_random_2<form_p_2]
            M2[new_bonds_2,9] = 1
            
            # -------------------------------------------------------------------------------------------------
            # dissociation of E-cad/Actin bond on the center region
            # -------------------------------------------------------------------------------------------------
            if (N_exist_bond_1!=0 or N_exist_bond_2!=0) and node_edge_yesorno[i] == 0:
                broken_p = 1.0/k5r_central*deltat*freq_update   # a constant
                broken_index_1 = exist_bond_1[np.random.rand(N_exist_bond_1)<broken_p] 
                broken_index_2 = exist_bond_2[np.random.rand(N_exist_bond_2)<broken_p]
                
                num_broken_bond_1[i] = broken_index_1.size
                num_broken_bond_2[i] = broken_index_2.size
                M1[broken_index_1,9] = 0  
                M2[broken_index_2,9] = 0  
                
        # overall forces on contact
        Force_central = 2*(num_myo_actin_1)*F_myo
        print('Force_central: {}, num_myo_actin: {}'.format(Force_central,num_myo_actin_1))
        
        # -------------------------------------------------------------------------------------------------
        # Calculate force from the edge & dissociation of E-cad-trans-actin 
        # -------------------------------------------------------------------------------------------------
        # adhesion energy from E-cad-trans dimers
        num_trans = np.count_nonzero(M1[:,2]+M2[:,2])
        F_adh = num_trans*g_trans*10**6*k*T      # N*E (J/m^2*e6)
       
        # calculate force toward the center of each node
        # update all E-cad-trans-actin dissociation on the edge
        # and only E-cad-actin bundled dissociation in the center
        # mode1 is the normal mode. Angle theta changes as a result of force balance
        # mode2: After the angle reaches a certain value, the angle will be fixed.

        edge_index = 0
        x_F_contact = Force_central-F_adh
        if x_F_contact > All_F_on_Ecad:
            memb_angle = 0
        else:
            memb_angle = np.arccos(x_F_contact/All_F_on_Ecad)

        for i in range(num_nodes_incenter):

            # index of lattice point in the ith element
            indexoflattice = np.array(matchedlist[i])   
            
            if node_edge_yesorno[i] == 1: 
                
                ####################################################################################
                # Calculation of the contact angle
                ####################################################################################
                # all exist E-cad-trans-actin on lattice field                    
                exist_bond_1 = indexoflattice[M1[indexoflattice,9]!=0]
                N_exist_bond_1 = exist_bond_1.size
                
                exist_bond_2 = indexoflattice[M2[indexoflattice,9]!=0]
                N_exist_bond_2 = exist_bond_2.size
                
                ####################################################################################
                # Calculate the E-cad actin dissociation probability
                ####################################################################################
                if (N_exist_bond_1 != 0) or (N_exist_bond_2 != 0):
                    
                    # count all lattice points on the nearby elements
                    all_lattice = lattice_nearby_node[edge_index]
                    
                    # Number of E_trans dimers on the nearby elements
                    N_Etrans_all_lattice = all_lattice[(M1[all_lattice,2]==1)].size
                    
                    # Force on the E-cad-trans dimer (lambda)
                    F_E_cad_trans = All_F_on_Ecad*np.sin(memb_angle)*(np.sum(Area_all[theta_nearby_node[edge_index]-1])/area_onedge)/N_Etrans_all_lattice
                    
                    # Force on the Factin = lambda/cos(theta)/sin(theta) (pN)
                    F_alpha_actin = F_E_cad_trans/(theta_cons)*10**-12
                    F_alpha_all = np.append(F_alpha_all,F_alpha_actin*10**12)
                    
                    if Force_sensitive == 1:    
                        broken_p = 1.0/alpha_dis(F_alpha_actin)*deltat*freq_update
                    elif Force_sensitive == 2:
                        broken_p = 1.0/k5r_central*deltat*freq_update
                    else:
                        print("Error in Force sensitive Mode")
                        break

                ####################################################################################
                # change E-cad/Actin unbinding rate
                ####################################################################################
                if (mode ==2) and (currentT > T_fix_angle) and (node_fixed_angle[edge_index] != 0):
                    k4_SF_temp = 0.02
                    # Update actin bundle dissociation rate
                    k4_1_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1,k4_SF_temp)
                    k4_2_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1,k4_SF_temp)
                    
                    k4_get = k4_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    k4_SF_tilt = np.append(k4_SF_tilt,k4_get)
                    
                edge_index += 1

                # -------------------------------------------------------------------------------------------------
                # dissociation of E-cad/Actin on the edge
                # -------------------------------------------------------------------------------------------------
                if N_exist_bond_1 !=0:
                    broken_index_1 = exist_bond_1[np.random.rand(N_exist_bond_1)<broken_p] 
                    
                    num_broken_bond_1[i] = broken_index_1.size
                    M1[broken_index_1,9] = 0                   # update M1 for broken E-trans-actin bonds
                
                if N_exist_bond_2 !=0:
                    broken_index_2 = exist_bond_2[np.random.rand(N_exist_bond_2)<broken_p] 
                    
                    num_broken_bond_2[i] = broken_index_2.size
                    M2[broken_index_2,9] = 0                   # update M1 for broken E-trans-actin bonds
            
            # -------------------------------------------------------------------------------------------------
            # Update E-cad concentration in the continuum model
            # -------------------------------------------------------------------------------------------------
            # Number of unbound E-trans dimers on the lattice field
            N_trans_unbound = indexoflattice[(M1[indexoflattice,2]==1)&(M1[indexoflattice,9]==0)&(M2[indexoflattice,9]==0)].size              # number of free trans dimers
            N_E_trans_actin_1[i] = indexoflattice[(M1[indexoflattice,2]==1)&(M1[indexoflattice,9]!=0)].size
            N_E_trans_actin_lattice_1 += N_E_trans_actin_1[i]
            
            N_E_trans_actin_2[i] = indexoflattice[(M2[indexoflattice,2]==1)&(M2[indexoflattice,9]!=0)].size
            N_E_trans_actin_lattice_2 += N_E_trans_actin_2[i]
            
            # Concentration of unbound E-trans dimers and E-trans-actin
            con_updated_unbound = N_trans_unbound/Area_all[i]
            con_bound_1 = N_E_trans_actin_1[i]/Area_all[i]
            con_bound_2 = N_E_trans_actin_2[i]/Area_all[i]
            
            # Update free E-trans dimer concentration
            E_trans_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1,con_updated_unbound)
            
            # Update E-cad-trans-actin concentration
            E_trans_actin_1_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1,con_bound_1)
            E_trans_actin_2_Field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1,con_bound_2)
        
        # -------------------------------------------------------------------------------------------------
        # Record concentration
        # -------------------------------------------------------------------------------------------------
        i_c = 0     # index of central nodes
        i_e = 0     # index of edge nodes
        
        
        if t%record_freq == 0: 
            con_E_trans_actin_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_E_trans_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_actin_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_actin_a_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_Rac_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_Rac_a_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_RhoA_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_RhoA_a_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_myo_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_myo_a_all = np.zeros(num_node_central+num_node_edge,dtype=float)
            con_myo_actin_all = np.zeros(num_node_central+num_node_edge,dtype=float)

            con_E_trans_actin_cent = np.zeros(num_node_central,dtype=float)
            con_E_trans_cent = np.zeros(num_node_central,dtype=float)
            con_actin_cent = np.zeros(num_node_central,dtype=float)
            con_actin_a_cent = np.zeros(num_node_central,dtype=float)
            con_Rac_cent = np.zeros(num_node_central,dtype=float)
            con_Rac_a_cent = np.zeros(num_node_central,dtype=float)
            con_RhoA_cent = np.zeros(num_node_central,dtype=float)
            con_RhoA_a_cent = np.zeros(num_node_central,dtype=float)
            con_myo_cent = np.zeros(num_node_central,dtype=float)
            con_myo_a_cent = np.zeros(num_node_central,dtype=float)
            con_myo_actin_cent = np.zeros(num_node_central,dtype=float)
            
            con_E_trans_actin_edge = np.zeros(num_node_edge,dtype=float)
            con_E_trans_edge = np.zeros(num_node_edge,dtype=float)
            con_actin_edge = np.zeros(num_node_edge,dtype=float)
            con_actin_a_edge = np.zeros(num_node_edge,dtype=float)
            con_Rac_edge = np.zeros(num_node_edge,dtype=float)
            con_Rac_a_edge = np.zeros(num_node_edge,dtype=float)
            con_RhoA_edge = np.zeros(num_node_edge,dtype=float)
            con_RhoA_a_edge = np.zeros(num_node_edge,dtype=float)
            con_myo_edge = np.zeros(num_node_edge,dtype=float)
            con_myo_a_edge = np.zeros(num_node_edge,dtype=float)
            con_myo_actin_edge = np.zeros(num_node_edge,dtype=float)

            con_Rac_2_cent = np.zeros(num_node_central,dtype=float)
            con_Rac_a_2_cent = np.zeros(num_node_central,dtype=float)
                    
            con_Rac_2_edge = np.zeros(num_node_edge,dtype=float)
            con_Rac_a_2_edge = np.zeros(num_node_edge,dtype=float)
                   
            np.zeros(num_node_edge,dtype=float)
            # Recording the number of mono, trans, cis on the lattice field 
            num_all_M1 = np.count_nonzero(M1[:,7])
            num_all_M2 = np.count_nonzero(M2[:,7])
            num_mono_M1 = np.count_nonzero(M1[:,1])
            num_mono_M2 = np.count_nonzero(M2[:,1])
            num_trans = np.count_nonzero(M1[:,2])
            num_trans_inring = np.count_nonzero(M1[inring_left,2])+np.count_nonzero(M1[inring_right,2])
            num_cis_incont_M1 = np.count_nonzero(M1[:,3]|M1[:,8])
            num_cis_incont_M2 = np.count_nonzero(M2[:,3]|M2[:,8])
            num_cad_actin_M1 = np.count_nonzero(M1[:,9])
            num_cad_actin_M2 = np.count_nonzero(M2[:,9])
            num_immobile_M1 = np.count_nonzero(M1[:,3]|M1[:,8]|M1[:,9])
            num_immobile_M2 = np.count_nonzero(M2[:,3]|M2[:,8]|M2[:,9])
            num_center_M1 = np.count_nonzero(M1[incent,7])
            num_center_M2 = np.count_nonzero(M2[incent,7])
            num_inringleft_M1 = np.count_nonzero(M1[inring_left,7])
            num_inringleft_M2 = np.count_nonzero(M2[inring_left,7])
            num_inringright_M1 = np.count_nonzero(M1[inring_right,7])
            num_inringright_M2 = np.count_nonzero(M2[inring_right,7])
            
            num_cad_all[t//record_freq,:] = (t//record_freq, num_all_M1, num_all_M2, num_mono_M1, num_mono_M2,
                                             num_trans, num_trans_inring, 
                                             num_cis_incont_M1, num_cis_incont_M2,  num_cad_actin_M1, num_cad_actin_M2,
                                             num_immobile_M1, num_immobile_M2, num_center_M1, num_center_M2,
                                             num_inringleft_M1, num_inringleft_M2, num_inringright_M1, num_inringright_M2)
            
            
            if mode == 2:
                if len(k10r_all) == 0:
                    k4_SF_tilt = k4_1
                    k4_SF_untilt = k4_1
                    con_myo_tilt = np.mean(con_myo_tilt)
                else:
                    k4_SF_tilt = k4_1
                    k4_SF_untilt = k4_1
                    con_myo_tilt = np.mean(con_myo_tilt)
            for i in range(num_nodes_incenter):
                if node_edge_yesorno[i] == 0: 
                    con_actin_cent[i_c] = actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_actin_a_cent[i_c] = actin_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_E_trans_actin_cent[i_c] = E_trans_actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_E_trans_cent[i_c] = E_trans_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_Rac_cent[i_c] = Rac_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_Rac_a_cent[i_c] = Rac_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_RhoA_cent[i_c] = RhoA_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_RhoA_a_cent[i_c] = RhoA_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_myo_cent[i_c] = myo_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_myo_a_cent[i_c] = myo_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_myo_actin_cent[i_c] = myo_actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    
                    con_Rac_2_cent[i_e] = Rac_2_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_Rac_a_2_cent[i_e] = Rac_a_2_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    
                    i_c += 1
                else:
                    con_actin_edge[i_e] = actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_actin_a_edge[i_e] = actin_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_E_trans_actin_edge[i_e] = E_trans_actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_E_trans_edge[i_e] = E_trans_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_Rac_edge[i_e] = Rac_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_Rac_a_edge[i_e] = Rac_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_RhoA_edge[i_e] = RhoA_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_RhoA_a_edge[i_e] = RhoA_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_myo_edge[i_e] = myo_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_myo_a_edge[i_e] = myo_a_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_myo_actin_edge[i_e] = myo_actin_1_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    
                    con_Rac_2_edge[i_e] = Rac_2_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    con_Rac_a_2_edge[i_e] = Rac_a_2_Field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,i+1,1)
                    
                    i_e += 1
                    
                    
                          
        #-------------------------------------------------------------------------------------------------
        # PDE SOLVE - TIME LOOP
        #-------------------------------------------------------------------------------------------------  

        rd_problem.Solve()
        print("solved")
        
        t_100 += 1         
        

        #-------------------------------------------------------------------------------------------------
        # CONTROL LOOP UPDATE - TIME LOOP
        #-------------------------------------------------------------------------------------------------        
        # start time, stop time, time increment 
        
        rd_controlLoop.TimesSet(currentT,currentT+deltat*freq_update,deltat*freq_update)
    
        #-------------------------------------------------------------------------------------------------
        # OUTPUT - TIME LOOP
        #-------------------------------------------------------------------------------------------------
        if t%10000 == 0:
             # Export results
            fields = iron.Fields()
            fields.CreateRegion(region)
        
            # Create new path
            newpath = r"./output"
            if not os.path.exists(newpath):
                os.makedirs(newpath)
    
            #Numbering the files from 0 to n_loops
            num_file = "{0:0" + str(len(str(t))) + "}"
        
    
            #Naming the file properly
            name = "output/AdherensJunctions_" + str(t)
    
            #Exporting the files
            fields.NodesExport(name,"FORTRAN")
            fields.ElementsExport(name,"FORTRAN")
            print("Writing", name, "of", str(t))
            
            fields.Finalise()
    # Record Concentrations
    if t%record_freq == 0:
        # Check for mass conservation
        con_actin_all = (sum(con_actin_cent*Area_all[node_central-1])+sum(con_actin_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_actin_a_temp_all = sum(con_actin_a_cent*Area_all[node_central-1])/sum(Area_all)
        con_E_trans_actin_all = (sum(con_E_trans_actin_cent*Area_all[node_central-1])+sum(con_E_trans_actin_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_E_trans_all = (sum(con_E_trans_cent*Area_all[node_central-1])+sum(con_E_trans_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_Rac_all = (sum(con_Rac_cent*Area_all[node_central-1])+sum(con_Rac_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_Rac_a_all = (sum(con_Rac_a_cent*Area_all[node_central-1])+sum(con_Rac_a_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_RhoA_all = (sum(con_RhoA_cent*Area_all[node_central-1])+sum(con_RhoA_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_RhoA_a_all = (sum(con_RhoA_a_cent*Area_all[node_central-1])+sum(con_RhoA_a_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_myo_all = (sum(con_myo_cent*Area_all[node_central-1])+sum(con_myo_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_myo_a_all = (sum(con_myo_a_cent*Area_all[node_central-1])+sum(con_myo_a_edge*Area_all[node_edge-1]))/sum(Area_all)
        con_myo_actin_all = sum(con_myo_actin_cent*Area_all[node_central-1])/sum(Area_all)
        
        # Average concentrations on the edge
        con_actin_cent_avg = sum(con_actin_cent)/num_node_central
        con_actin_a_cent_avg = sum(con_actin_a_cent)/num_node_central
        con_E_trans_actin_cent = sum(con_E_trans_actin_cent)/num_node_central
        con_E_trans_cent = sum(con_E_trans_cent)/num_node_central
        con_Rac_cent = sum(con_Rac_cent)/num_node_central
        con_Rac_a_cent = sum(con_Rac_a_cent)/num_node_central
        con_RhoA_cent = sum(con_RhoA_cent)/num_node_central
        con_RhoA_a_cent = sum(con_RhoA_a_cent)/num_node_central
        con_myo_cent = sum(con_myo_cent)/num_node_central
        con_myo_a_cent = sum(con_myo_a_cent)/num_node_central
        con_myo_actin_cent = sum(con_myo_actin_cent)/num_node_central
        
        con_Rac_2_cent = sum(con_Rac_2_cent)/num_node_central
        con_Rac_a_2_cent = sum(con_Rac_a_2_cent)/num_node_central

        # Average concentrations on the edge
        con_actin_edge_avg = sum(con_actin_edge)/num_node_edge
        con_actin_a_edge_avg = sum(con_actin_a_edge)/num_node_edge
        con_E_trans_actin_edge = sum(con_E_trans_actin_edge)/num_node_edge
        con_E_trans_edge = sum(con_E_trans_edge)/num_node_edge
        con_Rac_edge = sum(con_Rac_edge)/num_node_edge
        con_Rac_a_edge = sum(con_Rac_a_edge)/num_node_edge
        con_RhoA_edge = sum(con_RhoA_edge)/num_node_edge
        con_RhoA_a_edge = sum(con_RhoA_a_edge)/num_node_edge
        con_myo_edge = sum(con_myo_edge)/num_node_edge
        con_myo_a_edge = sum(con_myo_a_edge)/num_node_edge
        con_myo_actin_edge = sum(con_myo_actin_edge)/num_node_edge

        con_Rac_2_edge = sum(con_Rac_2_edge)/num_node_edge
        con_Rac_a_2_edge = sum(con_Rac_a_2_edge)/num_node_edge

        # Number of E-cad after
        num_bound_onedge_after = np.array(np.where((M1[onedge,2]==1)&(M1[onedge,9]!=0))).size
        num_bound_incenter_after = np.array(np.where((M1[incenter,2]==1)&(M1[incenter,9]!=0))).size
        
        # Concentration of E-cad
        con_bound_onedge = np.array(np.where((M1[onedge,2]==1)&(M1[onedge,9]!=0))).size/area_onedge
        con_bound_incenter = np.array(np.where((M1[incenter,2]==1)&(M1[incenter,9]!=0))).size/area_incenter

        #N_E_trans_actin_continuum = con_E_trans_actin*Area_all
        N_E_trans_actin_lattice = np.array(np.where((M1[:,2]==1)&(M1[:,9]!=0))).size
        
        file_con.write("#################################")
        file_con.write("Current Time: %0.4f\n\n"%currentT)
        #file_con.write("N_E_trans_actin_continuum: %0.4f\n"%sum(np.rint(N_E_trans_actin_continuum)))
        file_con.write("N_E_trans_actin_lattice: %0.4f\n\n"%N_E_trans_actin_lattice)

        file_con.write("ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL-ALL\n")
        
        file_con.write("Concentration of proteins ALL\n")
        file_con.write("Concentration of actin: %0.4f\n"%con_actin_all)
        file_con.write("Concentration of actin_a: %0.4f\n"%con_actin_a_temp_all)
        file_con.write("Concentration of E_trans_actin: %0.4f\n"%con_E_trans_actin_all)
        file_con.write("Concentration of E_trans_actin max: %0.4f\n"%np.amax(con_E_trans_actin_all))
        file_con.write("Concentration of E_trans: %0.4f\n"%con_E_trans_all)
        file_con.write("Concentration of E_trans max: %0.4f\n"%np.amax(con_E_trans_all))
        file_con.write("Concentration of Rac: %0.4f\n"%con_Rac_all)
        file_con.write("Concentration of Rac_a: %0.4f\n"%con_Rac_a_all)
        file_con.write("Concentration of RhoA: %0.4f\n"%con_RhoA_all)
        file_con.write("Concentration of RhoA_a: %0.4f\n"%con_RhoA_a_all)
        file_con.write("Concentration of myo: %0.4f\n"%con_myo_all)
        file_con.write("Concentration of myo_a: %0.4f\n"%con_myo_a_all)
        file_con.write("Concentration of myo_actin: %0.4f\n"%con_myo_actin_all)
        
        file_con.write("EDGE-EDGE-EDGE-EDGE-EDGE-EDGE-EDGE-EDGE-EDGE-EDGE-EDGE-EDGE\n")
        
        file_con.write("Concentration of proteins on edge\n")
        file_con.write("Concentration of actin: %0.4f\n"%con_actin_edge_avg)
        file_con.write("Concentration of actin_a: %0.4f\n"%con_actin_a_edge_avg)
        file_con.write("Concentration of E_trans_actin: %0.4f\n"%con_E_trans_actin_edge)
        file_con.write("Concentration of E_trans_actin max: %0.4f\n"%np.amax(con_E_trans_actin_edge))
        file_con.write("Concentration of E_trans: %0.4f\n"%con_E_trans_edge)
        file_con.write("Concentration of E_trans max: %0.4f\n"%np.amax(con_E_trans_edge))
        file_con.write("Concentration of Rac_1: %0.4f\n"%con_Rac_edge)
        file_con.write("Concentration of Rac_a_1: %0.4f\n"%con_Rac_a_edge)
        file_con.write("Concentration of Rac_2: %0.4f\n"%con_Rac_2_edge)
        file_con.write("Concentration of Rac_a_2: %0.4f\n"%con_Rac_a_2_edge)
        file_con.write("Concentration of RhoA: %0.4f\n"%con_RhoA_edge)
        file_con.write("Concentration of RhoA_a: %0.4f\n"%con_RhoA_a_edge)
        file_con.write("Concentration of myo: %0.4f\n"%con_myo_edge)
        file_con.write("Concentration of myo_a: %0.4f\n"%con_myo_a_edge)
        file_con.write("Concentration of myo_actin: %0.4f\n"%con_myo_actin_edge)
        
        file_con.write("CENT-CENT-CENT-CENT-CENT-CENT-CENT-CENT-CENT-CENT-CENT-CENT\n")
        
        file_con.write("Concentration of proteins in center\n")
        file_con.write("Concentration of actin: %0.4f\n"%con_actin_cent_avg)
        file_con.write("Concentration of actin_a: %0.4f\n"%con_actin_a_cent_avg)
        file_con.write("Concentration of E_trans_actin: %0.4f\n"%con_E_trans_actin_cent)
        file_con.write("Concentration of E_trans_actin max: %0.4f\n"%np.amax(con_E_trans_actin_cent))
        file_con.write("Concentration of E_trans: %0.4f\n"%con_E_trans_cent)
        file_con.write("Concentration of E_trans max: %0.4f\n"%np.amax(con_E_trans_cent))
        file_con.write("Concentration of Rac_1: %0.4f\n"%con_Rac_cent)
        file_con.write("Concentration of Rac_a_1: %0.4f\n"%con_Rac_a_cent)
        file_con.write("Concentration of Rac_2: %0.4f\n"%con_Rac_2_cent)
        file_con.write("Concentration of Rac_a_2: %0.4f\n"%con_Rac_a_2_cent)
        file_con.write("Concentration of RhoA: %0.4f\n"%con_RhoA_cent)
        file_con.write("Concentration of RhoA_a: %0.4f\n"%con_RhoA_a_cent)
        file_con.write("Concentration of myo: %0.4f\n"%con_myo_cent)
        file_con.write("Concentration of myo_a: %0.4f\n"%con_myo_a_cent)
        file_con.write("Concentration of myo_actin: %0.4f\n"%con_myo_actin_cent)
        
        #file_con.write("num_newbond: %0.4f\n"%sum(num_new_bond))     
        file_con.write("Average angle: %0.4f\n"%memb_angle)
        
        all_actin_right = np.mean(con_actin_edge[node_edge_right]+con_actin_a_edge[node_edge_right])
        all_actin_left = np.mean(con_actin_edge[node_edge_left]+con_actin_a_edge[node_edge_left])
        all_actin_cent = con_actin_cent_avg+con_actin_a_cent_avg
        file_con.write("Con E_trans_actin on left edge: %0.4f\n"%num_inringleft_M1)
        file_con.write("Con E_trans_actin on right edge: %0.4f\n"%num_inringright_M1)
        file_con.write("Con E_trans_actin in center: %0.4f\n"%num_center_M1)
        file_con.write("Con actin on left edge: %0.4f\n"%all_actin_left)
        file_con.write("Con actin on right edge: %0.4f\n"%all_actin_right)
        file_con.write("Con actin in center: %0.4f\n"%all_actin_cent)
        file_con.write("Force_central: {}, F_cort: {}".format(Force_central,F_cort))
        if F_alpha_all.size == 0:
            F_alpha_all = 0
        
        file_con.write("F_alpha_actin: %0.4f\n\n"%(np.mean(F_alpha_all)))
        
        edgeangle.append(memb_angle)
        
        
# Finalise OpenCMISS-Iron

np.savetxt('num_cad_all.txt',num_cad_all)
np.savetxt('edgeangle.txt',np.array(edgeangle))

file_M1.close()
file_M2.close()
file_con.close()
iron.Finalise()