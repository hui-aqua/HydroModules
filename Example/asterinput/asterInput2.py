loadr=[]               
loadr.append( _F(CHARGE=gF), )       
loadr.append( _F(CHARGE=fixed), )       
loadr.append( _F(CHARGE=sinkF1), )
        
if k == 1:
    resn = DYNA_NON_LINE(CARA_ELEM=elemprop,
                    CHAM_MATER=fieldmat,
                    COMPORTEMENT=(_F(DEFORMATION='GROT_GDEP',
                                     GROUP_MA=('twines','ropes'),
                                     RELATION='CABLE'),
                                  _F(DEFORMATION='PETIT',
                                     GROUP_MA=('topRings',),
                                     RELATION='ELAS'),
                                 ),
                    CONVERGENCE=_F(ITER_GLOB_MAXI=1000,
                                   RESI_GLOB_RELA=0.0001 ),
                    EXCIT=(loadr),
                    OBSERVATION=_F(GROUP_MA=('twines','ropes','topRings',),
                                    NOM_CHAM='DEPL',
                                    NOM_CMP=('DX','DY','DZ'),
                                    INST=k+dt,
                                    OBSE_ETAT_INIT='NON'),
                    SCHEMA_TEMPS=_F(FORMULATION='DEPLACEMENT',
                                   SCHEMA="HHT",
                                   ALPHA=-10.3,
                                   ),
                                   #add damping stablize the oscilations Need to study in the future
                    INCREMENT=_F(LIST_INST=times,INST_FIN=(k)*dt),
                    MODELE=model)
else:
    for i in range (NODEnumber):
        grpno = 'node%01g' %(i+1)
        l[i]=AFFE_CHAR_MECA( FORCE_NODALE=_F(GROUP_NO= (grpno),
                                             FX= Fnh[i][0],
                                             FY= Fnh[i][1],
                                             FZ= Fnh[i][2],
                                             ),
                             MODELE=model)
        
    for i in range (NODEnumber):
        loadr.append( _F(CHARGE=l[i],), )

    resn = DYNA_NON_LINE(CARA_ELEM=elemprop,
    				    CHAM_MATER=fieldmat,
    				    reuse=resn,
                        ETAT_INIT=_F(EVOL_NOLI=resn),
                        COMPORTEMENT=(_F(DEFORMATION='GROT_GDEP',
                                           GROUP_MA=('twines','ropes'),
                                           RELATION='CABLE'),
                                        _F(DEFORMATION='PETIT',
                                           GROUP_MA=('topRings',),
                                           RELATION='ELAS'),
                                       ),
                        CONVERGENCE=_F(ITER_GLOB_MAXI=1000,
                                       RESI_GLOB_RELA=0.0001 ),
                        EXCIT=(loadr),
                        OBSERVATION=_F(GROUP_MA=('twines','ropes','topRings',),
                                    NOM_CHAM='DEPL',
                                    NOM_CMP=('DX','DY','DZ'),
                                    INST=k+dt,
                                    OBSE_ETAT_INIT='NON'),
                        SCHEMA_TEMPS=_F(FORMULATION='DEPLACEMENT',
                                       SCHEMA="HHT",
                                        ALPHA=-10.3
                                       ),
                                       #add damping stablize the oscilations Need to study in the future
                        INCREMENT=_F(LIST_INST=times,INST_FIN=(k)*dt),
                        MODELE=model,
                        )
    if k < itimes:
        for i in range(NODEnumber):
            DETRUIRE(CONCEPT=_F(NOM=(l[i])))
            
tblp1 = POST_RELEVE_T(ACTION=(_F(OPERATION='EXTRACTION',      # For Extraction of values
                                 INTITULE='Nodal Displacements',    # Name of the table in .resu file
                                 # The result from which values will be extracted(STAT_NON_LINE)
                                 RESULTAT=resn,
                                 # Field to extract. DEPL = Displacements
                                 NOM_CHAM=('DEPL'),
                                 # TOUT_CMP='OUI',
                                 # Components of DISP to extract
                                 NOM_CMP=('DX', 'DY', 'DZ'),
                                 GROUP_NO='allnodes',               # Extract only for nodes of group DISP
                                 # STAT_NON_LINE calculates for 10 INST. I want only last INST
                                 INST=(k)*dt,
                                 ),),
                      )
tblp2 = POST_RELEVE_T(ACTION=(_F(OPERATION='EXTRACTION',      # For Extraction of values
                                 INTITULE='Nodal Displacements',    # Name of the table in .resu file
                                 # The result from which values will be extracted(STAT_NON_LINE)
                                 RESULTAT=resn,
                                 # Field to extract. VITE = velocity,
                                 NOM_CHAM=('VITE'),
                                 # TOUT_CMP='OUI',
                                 # Components of DISP to extract
                                 NOM_CMP=('DX', 'DY', 'DZ'),
                                 GROUP_NO='allnodes',               # Extract only for nodes of group DISP
                                 # STAT_NON_LINE calculates for 10 INST. I want only last INST
                                 INST=(k)*dt,
                                 ),),
                      )

tblp3 = POST_RELEVE_T(ACTION=(_F(OPERATION='EXTRACTION',      # For Extraction of values
                                 INTITULE='Nodal Displacements',    # Name of the table in .resu file
                                 # The result from which values will be extracted(STAT_NON_LINE)
                                 RESULTAT=resn,
                                 # Field to extract. VITE = velocity,
                                 NOM_CHAM=('ACCE'),
                                 # TOUT_CMP='OUI',
                                 # Components of DISP to extract
                                 NOM_CMP=('DX', 'DY', 'DZ'),
                                 GROUP_NO='allnodes',               # Extract only for nodes of group DISP
                                 # STAT_NON_LINE calculates for 10 INST. I want only last INST
                                 INST=(k)*dt,
                                 ),),
                      )

posi_nodes = hdm.aster.get_position_aster(tblp1)
velo_nodes = hdm.aster.get_velocity_aster(tblp2)
acce_nodes = hdm.aster.get_acceleration_aster(tblp3)

DETRUIRE(CONCEPT=_F(NOM=(tblp1)))
DETRUIRE(CONCEPT=_F(NOM=(tblp2)))
DETRUIRE(CONCEPT=_F(NOM=(tblp3)))

timeFE = dt*k

u_current = min(1.0,(timeFE)/10.0)*np.ones((len(netting.output_hydro_element()), 3))*Uinput*wake.reduction_factors
# print(u_current)

load_netting[str(round(timeFE,3))]=(netting.cal_buoy_force(posi_nodes) + \
                                    netting.force_on_element(posi_nodes, u_current)).tolist()
load_pipes_top[str(round(timeFE,3))]=(pipes_top.cal_buoy_force(posi_nodes) + \
                                      pipes_top.force_on_element(posi_nodes, Uinput)).tolist()
load_rope[str(round(timeFE,3))]=(rope.cal_buoy_force(posi_nodes) + \
                                 rope.force_on_element(posi_nodes, Uinput)).tolist()

Fnh = netting.distribute_force(posi_nodes)+ \
      pipes_top.distribute_force(posi_nodes)+ \
      rope.distribute_force(posi_nodes)

Nodes_loads[round(timeFE,3)]=Fnh.tolist()
Nodes_velocity[round(timeFE,3)]=velo_nodes.tolist()
Nodes_position[round(timeFE,3)]=posi_nodes.tolist()


###>>> for easy post process
total_results['time'].append(round(timeFE,3))

total_results['total_force_on_net'].append((netting.distribute_force(posi_nodes).sum(axis=0)).tolist())
total_results['dynamic_force_on_net'].append((netting.force_on_element(posi_nodes, u_current).sum(axis=0)).tolist())
total_results['static_force_on_net'].append((netting.cal_buoy_force(posi_nodes).sum(axis=0)).tolist())

total_results['total_force_on_cable'].append((rope.distribute_force(posi_nodes).sum(axis=0)).tolist())
total_results['dynamic_force_on_cable'].append((rope.force_on_element(posi_nodes,Uinput).sum(axis=0)).tolist())
total_results['static_force_on_cable'].append((rope.cal_buoy_force(posi_nodes).sum(axis=0)).tolist())


total_results['total_force_on_top_pipe'].append((pipes_top.distribute_force(posi_nodes).sum(axis=0)).tolist())
total_results['dynamic_force_on_top_pipe'].append((pipes_top.force_on_element(posi_nodes,Uinput).sum(axis=0)).tolist())
total_results['static_force_on_top_pipe'].append((pipes_top.cal_buoy_force(posi_nodes).sum(axis=0)).tolist())

total_results['total_force'].append((Fnh.sum(axis=0)).tolist())

total_results['volume'].append(wake.cal_cage_volume(posi_nodes))

total_results['position_range_x'].append([np.min(posi_nodes[:,0]),np.max(posi_nodes[:,0])])
total_results['position_range_y'].append([np.min(posi_nodes[:,1]),np.max(posi_nodes[:,1])])
total_results['position_range_z'].append([np.min(posi_nodes[:,2]),np.max(posi_nodes[:,2])])

total_results['net_position_range_x'].append([np.min(posi_nodes[node_net][:,0]),np.max(posi_nodes[node_net][:,0])])
total_results['net_position_range_y'].append([np.min(posi_nodes[node_net][:,1]),np.max(posi_nodes[node_net][:,1])])
total_results['net_position_range_z'].append([np.min(posi_nodes[node_net][:,2]),np.max(posi_nodes[node_net][:,2])])



## >>> write results
with open(os.path.join(cwd,'pythonOutput','total_results.json'), "w") as json_file:
   json.dump(total_results,json_file)
json_file.close()               


with open(os.path.join(cwd,'pythonOutput','load_netting.json'), "w") as json_file:
   json.dump(load_netting,json_file)
json_file.close()
with open(os.path.join(cwd,'pythonOutput','load_pipes_top.json'), "w") as json_file:
   json.dump(load_pipes_top,json_file)
json_file.close()
with open(os.path.join(cwd,'pythonOutput','load_rope.json'), "w") as json_file:
   json.dump(load_rope,json_file)
json_file.close()

with open(os.path.join(cwd,'pythonOutput','Nodes_loads.json'), "w") as json_file:
   json.dump(Nodes_loads,json_file)
json_file.close()
with open(os.path.join(cwd,'pythonOutput','Nodes_velocity.json'), "w") as json_file:
   json.dump(Nodes_velocity,json_file)
json_file.close()
with open(os.path.join(cwd,'pythonOutput','Nodes_position.json'), "w") as json_file:
   json.dump(Nodes_position,json_file)
json_file.close()
