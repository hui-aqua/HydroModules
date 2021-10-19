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
import json
import numpy as np
from numpy import pi
point=[[0,0,0],
       [1,1,1],
       [0,0,3]]

con_all=[[1,2],
         [2,3]]
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
    nodeID = Mesh_1.AddNode(float(each_node[0]), float(
        each_node[1]), float(each_node[2]))

for each_con in con_all:
    edge = Mesh_1.AddEdge([each_con[0], each_con[1]])

isDone = Mesh_1.Compute()

# naming  the group
# naming the node
# GROUP_NO
allnodes = Mesh_1.CreateEmptyGroup(SMESH.NODE, 'allnodes')
nbAdd = allnodes.AddFrom(Mesh_1.GetMesh())
smesh.SetName(allnodes, 'allnodes')

meshname='test.med'

try:
    Mesh_1.ExportMED(os.path.join(os.getcwd(), meshname))
    pass
except:
    print('ExportMED() failed. Invalid file name?')

killMyPort(os.getenv('NSPORT'))
meshinfo={
    "node":point,
    "elem":con_all,
}
with open(os.path.join(os.getcwd(), 'meshInformation.json'), 'w') as json_file:
    json.dump(meshinfo, json_file,indent=4)
json_file.close()

