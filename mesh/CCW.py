'''
/--------------------------------\
|    University of Stavanger     |
|           Hui Cheng            |
\--------------------------------/
Any questions about this code, please email: hui.cheng@uis.no
The center of the floating collar is (0,0,0)
Fish cage is along the Z- direction
Z=0 is the free surface
Z<0 is the water zone
Z>0 is the air zone
The sinkers are attached to the floating collar
'''

from killSalomeWithPort import killMyPort
from salome.smesh import smeshBuilder
import SMESH
import salome

import os
import sys
import ast
import json
import numpy as np
from numpy import pi


# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter template start
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# default value
parameters={
'Environment':
{
    'current':[0.5,0,0],   # current velcity [m/s]
    'fluidDensity':1025,  # density of fresh water/sea water [kg/m3]
    'waterDepth':60,
},
'CageShape':
  {
    'origin':[[0,0,0]],
    'cageDiameter':50, # [m] diameter of the fish cage
    'cageHeight':15,  # [m] height of the fish cage
    'cageConeHeight':3,
    'elementOverCir':36,
    'elementOverHeight':6,
    'elementOverCone':5,
  },
'Frame':
    {
    'topringDiameter':51,  #[m]
    'topringDensity':110,
    'bottomringDiameter':51,
    'bottomringDepth':16,
    'bottomringDensity':2400,
    'elementOverRings':36, # this value is better to be same with 'elementOverCir'.
    'elementOverHeight': 8,
    'topring_sec_dia':0.35,
    'topring_thickness':0.0185,
    'bottomring_sec_dia':0.25,
    'bottomring_thickness':0.0185,
    'rope_sec_dia':0.05,
    'weight_tip':980,
    },    
'Net':
  {
    'nettingType':'square',
    'Sn': 0.25,   # solidity ratio
    'twineDiameter': 1.5e-3, # [m]the twine diameter of the physical net
    'meshLength': 12e-3, # [m]the half mesh length
    'netYoungmodule':2e9, # [Pa]
    'netRho':1140.0, #[kg/m3] density of the net material
  },
}


# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter template finish
# #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

path_to_setting=os.path.join(os.getcwd(),'setting.py')
if os.path.isfile(path_to_setting):
    print('\nYes, the setting file exists. The default parameters will be overwritted. \n')
    with open(path_to_setting, 'r') as f:
        content = f.read()
        parameters = ast.literal_eval(content)
else:
    print("Using the default parameters for modelling!!\n")

floater_center = parameters['CageShape']['origin']
cr_top = float(parameters['CageShape']['cageCircumference']/pi/2.0)
cr_bottom = cr_top*float(parameters['CageShape']['bottom_top_ratio'])
dr=cr_top-cr_bottom  # 
cage_height = float(parameters['CageShape']['cageHeight'])
cage_cone_height = float(parameters['CageShape']['cageConeHeight'])

NT = int(parameters['CageShape']['cageCircumference']/parameters['CageShape']['element_length'])     # Number of the nodes in circumference
NN = int(cage_height/parameters['CageShape']['element_length'])  # number of section in the height, thus, the nodes should be NN+1
BN = int(cr_bottom/parameters['CageShape']['element_length']/np.sqrt(2))    # number of section along the cone, thus, the nodes should be NN+1

# below is the final information for outputing
# con_pipe 
# con_rope 
# con_netting 
# sur_netting 
# point


# generate the point coordinates matrix for cylindrical netting
# procedure for generating point_netting:
# 1. Generate the node on one cage based on rotating the nodes on one side using transcope matrix.
# 2. Generate the node on (all) netting through translation.
point_one_cage=[]
# Step 1. 
for i in range(0, NT):
    for j in range(NN):
        point_one_cage.append(
            [(cr_top-dr*j/float(NN)) * np.cos(i * 2 * pi / float(NT)),
             (cr_top-dr*j/float(NN)) * np.sin(i * 2 * pi / float(NT)),
              - j * cage_height / float(NN)])
    for j in range(BN):
        point_one_cage.append(
            [cr_bottom * ((BN - j) / BN) * np.cos(i * 2 * pi / float(NT)),
             cr_bottom * ((BN - j) / BN) * np.sin(i * 2 * pi / float(NT)),
             - cage_height - j * (cage_cone_height) / float(BN)])
point_one_cage.append([0, 0, -cage_cone_height-cage_height])  # the last point should be at the cone tip
      
number_of_point_one_netting= len(point_one_cage)
print(number_of_point_one_netting)    

# generate the point coordinates matrix for mooring part
r_top=cr_top+1
NN_frame=NN
BN_frame=NN
NT_frame=NT
r_bottom=cr_bottom+1

for i in range(0, NT_frame):
    for j in range(0, NN_frame):
        point_one_cage.append(
            [(r_top-dr*j/float(NN_frame)) * np.cos(i * 2 * pi / float(NT_frame)),
             (r_top-dr*j/float(NN_frame)) * np.sin(i * 2 * pi / float(NT_frame)),
             -cage_height*j/float(NN_frame)])
    for j in range(BN_frame):
        point_one_cage.append(
            [r_bottom * ((BN_frame - j) / BN_frame) * np.cos(i * 2 * pi / float(NT)),
             r_bottom * ((BN_frame - j) / BN_frame) * np.sin(i * 2 * pi / float(NT)),
             - cage_height - j * (float(parameters['Frame']['rope_depth_height_ratio'])*cage_height) / float(BN_frame)])
point_one_cage.append([0, 0, -float(parameters['Frame']['rope_depth_height_ratio'])*cage_height-cage_height])  # the last point should be at the cone tip             

number_of_point_one_frame=len(point_one_cage)-number_of_point_one_netting
number_of_point_one_cage=len(point_one_cage)
print(number_of_point_one_cage)
# Step 2        
point=[]
for center in floater_center:
    point+=(np.array(point_one_cage) + np.ones((len(point_one_cage), 3))*np.array(center)).tolist()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Generate connection and sur
# similar with the nodes, first is for one side, then for one cage, finally for all the cages.
# 1.one cage

con_in_one_netting = []
sur_in_one_netting = []
sur_in_one_netting_tri=[]
# horizontal con for netting
for i in range(BN+NN):
    con_in_one_netting.append([1+i,1+i+(BN+NN)*(NT-1)])
    for j in range(NT-1):
        con_in_one_netting.append([1+i+j*(BN+NN),1+i+(j+1)*(BN+NN)])

# vertical con for netting
for i in range(BN+NN-1):
    for j in range(NT):
        con_in_one_netting.append([1+i+j*(BN+NN),1+1+i+j*(BN+NN)])
        if j < NT-1:
            sur_in_one_netting.append([1+i+j*(BN+NN),1+1+i+j*(BN+NN),
                                    1+i+(j+1)*(BN+NN),1+1+i+(j+1)*(BN+NN)])
        else:
            pass
            sur_in_one_netting.append([1+i+j*(BN+NN),1+1+i+j*(BN+NN),
                                    1+i,1+1+i])                
# print(sur_in_one_cage)                #test    
# bottom part, triangular element
for j in range(NT):
    con_in_one_netting.append([number_of_point_one_netting,
                            (j+1)*(BN+NN)])
    if j < NT-1:
        sur_in_one_netting_tri.append([number_of_point_one_netting,
                                    (j+1)*(BN+NN),
                                    (j+1+1)*(BN+NN)
                                    ])
    else:
        sur_in_one_netting_tri.append([number_of_point_one_netting,
                                    BN+NN,
                                    NT*(BN+NN)
                                    ])

# generate the connection for frame
# print(len(point_netting))             #test
# print(number_of_point_one_frame)           #test
# con_in_one_frame = []
con_in_one_frame_top=[]
con_in_one_frame_line=[]

for j in range(NT_frame):
    # top
    con_in_one_frame_top.append([1+j*(BN+NN),1+j*(NN_frame+BN_frame)+number_of_point_one_netting])
    # circular
    if j<NT_frame-1:
        con_in_one_frame_top.append([1+j*(NN_frame+BN_frame)+number_of_point_one_netting,1+(j+1)*(NN_frame+BN_frame)+number_of_point_one_netting])
    else:
        con_in_one_frame_top.append([1+j*(NN_frame+BN_frame)+number_of_point_one_netting,1+number_of_point_one_netting])            
    #vertical lines
    for i in range(NN_frame+BN_frame-1):
        con_in_one_frame_line.append([1+i+j*(NN_frame+BN_frame)+number_of_point_one_netting,1+1+i+j*(NN_frame+BN_frame)+number_of_point_one_netting])
    con_in_one_frame_line.append([NN_frame+BN_frame+j*(NN_frame+BN_frame)+number_of_point_one_netting,number_of_point_one_cage])
    # TODO 
    for num in range(int(NN/2)):
        con_in_one_frame_line.append([1+NN+j*(BN+NN)-num*2,1+NN_frame+j*(NN_frame+BN_frame)+number_of_point_one_netting-2*num])

# Step 2  
sur_netting_tri=[]
con_one_cage=con_in_one_netting+con_in_one_frame_line+con_in_one_frame_top
num_of_con_in_one_cage=len(con_one_cage)


con_all=[]
con_pipe_top = []
con_rope = []
con_netting = []  # lines on netting
sur_netting=[]

for num in range(len(floater_center)):
    sur_netting     += (np.array(sur_in_one_netting)      + num*number_of_point_one_cage).tolist()
    sur_netting_tri += (np.array(sur_in_one_netting_tri)  + num*number_of_point_one_cage).tolist()
    con_all         += (np.array(con_one_cage)            + num*number_of_point_one_cage).tolist()
    con_pipe_top    += (np.array(con_in_one_frame_top)    + num*number_of_point_one_cage).tolist()
    con_rope        += (np.array(con_in_one_frame_line)   + num*number_of_point_one_cage).tolist()
    con_netting     += (np.array(con_in_one_netting)      + num*number_of_point_one_cage).tolist()
sur_netting +=  sur_netting_tri




# ###############salome########################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the below is the command in the Mesh, Salome.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ###############salome########################

salome.salome_init()
theStudy = salome.myStudy

smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh()

# add the pints into geometry
for each_node in point:
    nodeID = Mesh_1.AddNode(float(each_node[0]), float(each_node[1]), float(each_node[2]))

for each_con in con_all:
    edge = Mesh_1.AddEdge([each_con[0], each_con[1]])

isDone = Mesh_1.Compute()

# naming  the group
# naming the node
# GROUP_NO
allnodes = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'allnodes')
nbAdd = allnodes.AddFrom(Mesh_1.GetMesh())
smesh.SetName(allnodes, 'allnodes')

# generate the name for each node to assign the hydrodynamic forces.
for i in range(1, len(point) + 1):
    node1 = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'node%s' % i)
    nbAdd = node1.Add([i])
    smesh.SetName(node1, 'node{}'.format(str(i)))

# define the nodes on net top
net_top = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'net_top')
for num in range(len(floater_center)):
    index_start=number_of_point_one_cage*num
    nbAdd = net_top.Add([i for i in range(index_start+1, index_start+number_of_point_one_netting,BN+NN)])
smesh.SetName(net_top, 'net_top')

# define the nodes on  net bottom
net_bottom = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'net_bottom')
for num in range(len(floater_center)):
    index_start=number_of_point_one_cage*num+NN
    nbAdd = net_bottom.Add([i for i in range(index_start+1, index_start+number_of_point_one_netting,BN+NN)])
smesh.SetName(net_bottom, 'net_bottom')

# define the nodes on top ring
ring_top = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'ring_top')
for num in range(len(floater_center)):
    index_start=number_of_point_one_netting+number_of_point_one_cage*num
    nbAdd = ring_top.Add([i for i in range(index_start+1, index_start+number_of_point_one_frame,NN_frame+BN_frame)])
smesh.SetName(ring_top, 'ring_top')

# define the nodes on bottom ring
# ring_bottom = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'ring_bottom')
# for num in range(len(floater_center)):
#     index_start=number_of_point_one_netting+number_of_point_one_cage*num+NN_frame
#     nbAdd = ring_bottom.Add([i for i in range(index_start+1, index_start+number_of_point_one_frame+1,NN_frame+1)])
# smesh.SetName(ring_bottom, 'ring_bottom')

# define the bottom tips
bottom_tip = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'bottom_tip')
for num in range(len(floater_center)+1):
    nbAdd = bottom_tip.Add([number_of_point_one_netting+number_of_point_one_cage*num])
smesh.SetName(bottom_tip, 'bottom_tip')
 
# define the bottom tips
weight_tip = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'weight_tip')
for num in range(len(floater_center)+1):
    nbAdd = weight_tip.Add([number_of_point_one_cage+number_of_point_one_cage*num])
smesh.SetName(weight_tip, 'weight_tip')
 
# GROUP_MA
# defaults names for all the twines.
twines = Mesh_1.CreateEmptyGroup(SMESH.EDGE, 'twines')
for num in range(len(floater_center)):       
    index_start=num_of_con_in_one_cage*num
    nbAdd = twines.Add([i for i in range(index_start+1+NT, index_start+1+len(con_in_one_netting))])
smesh.SetName(twines, 'twines')

# defaults names for rope
ropes = Mesh_1.CreateEmptyGroup(SMESH.EDGE, 'ropes')
for num in range(len(floater_center)):       
    index_start=len(con_in_one_netting)+num_of_con_in_one_cage*num
    nbAdd = ropes.Add([i for i in range(index_start+1, index_start+1+len(con_in_one_frame_line))])
smesh.SetName(ropes, 'ropes')


# below should be set as HDPE beam
# # defaults names for top ring
topRings = Mesh_1.CreateEmptyGroup(SMESH.EDGE, 'topRings')
for num in range(len(floater_center)):
    index_start=len(con_in_one_netting)+len(con_in_one_frame_line)+num_of_con_in_one_cage*num
    nbAdd = topRings.Add([i for i in range(index_start+1, index_start+1+len(con_in_one_frame_top))])
    nbAdd = topRings.Add([i for i in range(num*num_of_con_in_one_cage+1, num*num_of_con_in_one_cage+1+NT)])
smesh.SetName(topRings, 'topRings')

# # defaults names for bottom ring
# bottomRings = Mesh_1.CreateEmptyGroup(SMESH.EDGE, 'bottomRings')
# for num in range(len(floater_center)):
#     index_start=len(con_in_one_netting)+len(con_in_one_frame_line)+len(con_in_one_frame_top)+num_of_con_in_one_cage*num
#     nbAdd = bottomRings.Add([i for i in range(index_start+1, index_start+1+len(con_in_one_frame_bottom))])
# smesh.SetName(bottomRings, 'bottomRings')


print('\n\nCase setup helper>>>>>>>>>>>>')
total_length_of_twine=float(smesh.GetLength(twines))
number_of_twine= len(con_in_one_netting*len(floater_center))
Ln_average=total_length_of_twine/number_of_twine
lp=float(parameters['Net']['meshLength'])
dws=float(parameters['Net']['twineDiameter'])*np.sqrt(Ln_average/lp)
dwh=float(parameters['Net']['twineDiameter'])*(Ln_average/lp)
section_twine = 0.25*np.pi*pow(dws, 2)
total_volume_of_twine = float(section_twine*total_length_of_twine)

### write to the meshinfo:
fb_netting=total_volume_of_twine*9.81*parameters['Environment']['fluidDensity']/number_of_point_one_netting




# give a name to the mesh
meshname = 'CCW.med'

meshinfo = {
    'Lines_rope': con_rope,
    'numberOfLines_rope': len(con_rope),
    'Lines_pipe_top': con_pipe_top,
    'numberOfLines_pipe_top': len(con_pipe_top),
    'Lines_netting': con_netting,
    'surfs_netting': sur_netting,
    'numberOfLines_netting': len(con_netting),
    'numberOfsurfs_netting': len(sur_netting),
    'Nodes':point,
    'numberOfNode':len(point),
    'meshName': meshname,
    'fb_netting':fb_netting,
    'dws':dws,
    'dwh':dwh,
    'caseInfo':parameters,
    'NT':NT,
    'NN':NN,
    'BN':BN,
}


try:
    Mesh_1.ExportMED(os.path.join(os.getcwd(), meshname))
    pass
except:
    print('ExportMED() failed. Invalid file name?')

killMyPort(os.getenv('NSPORT'))
print('\n'
      '  --------------------------------------\n'
      '  --     University of Stavanger      --\n'
      '  --         Hui Cheng (PhD student)  --\n'
      '  --       Lin Li (Medveileder)       --\n'
      '  --  Prof. Muk Chen Ong (supervisor) --\n'
      '  --------------------------------------\n'
      '  Any questions about this code,\n'
      '  please email: hui.cheng@uis.no\n')

print('<<<<<<<<<< Mesh Information >>>>>>>>>>')
print('Number of node is ' + str(len(point)) + '.')
print('Number of line is ' + str(len(con_all)) + '.')
print('Mesh file is generated and stored at '+str(os.path.join(os.getcwd(),meshname)))



with open(os.path.join(os.getcwd(), 'meshInformation.json'), 'w') as json_file:
    json.dump(meshinfo, json_file,indent=4)
json_file.close()


if __name__ == '__main__':
    pass
