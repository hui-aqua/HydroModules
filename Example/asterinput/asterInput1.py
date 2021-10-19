# ----------------------------------
# --   University of Stavanger    --
# --           Hui Cheng          --
# ----------------------------------
# Any questions about this code,
# please email: hui.cheng@uis.no
import sys
import os
import numpy as np
import json
import itertools

cwd = "PATH_OF_THE_CURRENT_FOLDER"
sys.path.append(os.path.join(cwd, 'asterinput/module'))
import hydroModules as hdm
import enviromentModules as env

# read meshinfo
with open(os.path.join(cwd, 'meshInformation.json')) as json_file:
    meshinfo = json.load(json_file)
json_file.close()

with open(os.path.join(cwd, 'setting.json')) as json_file:
    caseinfo = json.load(json_file)
json_file.close()


dt = 0.1
duration = 6
itimes = int(duration/dt)
tend = itimes*dt
Uinput = caseinfo['Environment']['current']
NODEnumber = meshinfo['numberOfNode']
Fnh = [[0.0, 0.0, 0.0]]*NODEnumber
l = ['None']*(NODEnumber)


# define hydrodynamic objects
netting = hdm.two_dimensional.netting("S3",
                                      meshinfo['surfs_netting'],
                                      caseinfo['Net']['Sn'],
                                      caseinfo['Net']['twineDiameter'],
                                      )
pipes_top = hdm.one_dimensional.pipe("M4",
                                     meshinfo['Lines_pipe_top'],
                                     caseinfo['Frame']['topring_sec_dia'],
                                     caseinfo['Frame']['topring_thickness'],
                                     permeability=False
                                     )
rope = hdm.one_dimensional.line("M4",
                                meshinfo['Lines_rope'],
                                0,
                                caseinfo['Frame']['rope_sec_dia'],
                                caseinfo['Frame']['rope_sec_dia'])

wake = hdm.weak_effect.net2netWeak('factor-0.8',
                                   meshinfo['Nodes'],
                                   netting.output_hydro_element(),
                                   Uinput,
                                   caseinfo['Net']['Sn'])

elements = {}
elements['netting'] = netting.output_hydro_element()
elements['pipes_top'] = pipes_top.output_hydro_element()
elements['rope'] = rope.output_hydro_element()

flat_list = itertools.chain(*elements['netting'])
node_net=list(set(flat_list))

# output info for the hydrodynamic elements
with open(os.path.join(cwd, 'asteroutput', 'hydro_elements.json'), "w") as json_file:
    json.dump(elements, json_file, indent=4)
json_file.close()



total_results = {'total_weight': caseinfo['Frame']['weight_per_metter']*\
                                 caseinfo['CageShape']['cageCircumference']*9.81,
                 'total_net_area':wake.cal_cage_net_area(np.array(meshinfo['Nodes'])),
                 
                 'time': [0.0],
                 
                 'total_force_on_net':   [[0.0, 0.0, 0.0]],
                 'dynamic_force_on_net': [[0.0, 0.0, 0.0]],
                 'static_force_on_net':  [[0.0, 0.0, 0.0]],
                 
                 'total_force_on_cable':   [[0.0, 0.0, 0.0]],
                 'dynamic_force_on_cable': [[0.0, 0.0, 0.0]],
                 'static_force_on_cable':  [[0.0, 0.0, 0.0]],
                 
                 'total_force_on_top_pipe': [[0.0, 0.0, 0.0]],
                 'dynamic_force_on_top_pipe': [[0.0, 0.0, 0.0]],
                 'static_force_on_top_pipe': [[0.0, 0.0, 0.0]],
                                 
                 'total_force': [[0.0, 0.0, 0.0]],
                 
                 'volume': [wake.cal_cage_volume(meshinfo['Nodes'])],
                 
                 'position_range_x': [[np.min(np.array(meshinfo['Nodes'])[:, 0]),
                                       np.max(np.array(meshinfo['Nodes'])[:, 0])]],
                 'position_range_y': [[np.min(np.array(meshinfo['Nodes'])[:, 1]),
                                       np.max(np.array(meshinfo['Nodes'])[:, 1])]],
                 'position_range_z': [[np.min(np.array(meshinfo['Nodes'])[:, 2]),
                                       np.max(np.array(meshinfo['Nodes'])[:, 2])]],
                 
                 'net_position_range_x': [[np.min(np.array(meshinfo['Nodes'])[node_net][:, 0]),
                                           np.max(np.array(meshinfo['Nodes'])[node_net][:, 0])]],
                 'net_position_range_y': [[np.min(np.array(meshinfo['Nodes'])[node_net][:, 1]),
                                           np.max(np.array(meshinfo['Nodes'])[node_net][:, 1])]],
                 'net_position_range_z': [[np.min(np.array(meshinfo['Nodes'])[node_net][:, 2]),
                                           np.max(np.array(meshinfo['Nodes'])[node_net][:, 2])]],
                 }


load_netting = {}
load_pipes_top = {}
load_rope = {}
Nodes_loads = {}
Nodes_velocity = {}
Nodes_position = {'initial': meshinfo['Nodes']}


DEBUT(PAR_LOT='NON',
      IGNORE_ALARM=("SUPERVIS_25", "DISCRETE_26", "UTILITAI8_56")
      )

# read mesh
mesh = LIRE_MAILLAGE(UNITE=20)

# define element type
model = AFFE_MODELE(AFFE=(_F(GROUP_MA=('twines', 'ropes'),
                             MODELISATION=('CABLE'),
                             PHENOMENE='MECANIQUE'),
                          _F(GROUP_MA=('topRings',),
                             MODELISATION=('POU_D_E'),
                             PHENOMENE='MECANIQUE')
                          ),
                    MAILLAGE=mesh)

# define element geometrical properties
elemprop = AFFE_CARA_ELEM(CABLE=(_F(GROUP_MA=('twines'),
                                    N_INIT=10.0,
                                    SECTION=0.25*np.pi*pow(meshinfo['dws'], 2)),  # 10twines as a seam
                                 _F(GROUP_MA=('ropes'),
                                    N_INIT=40.0,
                                    SECTION=0.25*np.pi*pow(caseinfo['Frame']['rope_sec_dia'], 2)),
                                 ),
                          POUTRE=(_F(GROUP_MA=('topRings'),
                                     SECTION='CERCLE',
                                     CARA=('R', 'EP'),
                                     VALE=(caseinfo['Frame']['topring_sec_dia']/2,
                                           caseinfo['Frame']['topring_thickness'])),
                                  ),
                          MODELE=model)
# material base
# netting according to SCTsetting.py
net = DEFI_MATERIAU(CABLE=_F(EC_SUR_E=0.0001),
                    ELAS=_F(E=caseinfo['Net']['netYoungmodule'],
                            NU=0.3,  # No meaning for 1D element
                            RHO=caseinfo['Net']['netRho']))
# PE rope
pe = DEFI_MATERIAU(CABLE=_F(EC_SUR_E=0.0001),
                   ELAS=_F(E=1e9,
                           NU=0.2,
                           RHO=1100.0))
# hdpe
hdpe = DEFI_MATERIAU(ELAS=_F(E=3e9,
                             NU=0.2,
                             RHO=958.0))

fieldmat = AFFE_MATERIAU(AFFE=(_F(GROUP_MA=('twines',),
                                  MATER=(net)),
                               _F(GROUP_MA=('ropes',),
                                  MATER=(pe)),
                               _F(GROUP_MA=('topRings'),
                                  MATER=(hdpe)),
                               ),
                         MODELE=model)

# load
gF = AFFE_CHAR_MECA(PESANTEUR=_F(DIRECTION=(0.0, 0.0, -1.0),
                                 GRAVITE=9.81,
                                 GROUP_MA=('twines', 'ropes', 'topRings')),
                    MODELE=model)

fixed = AFFE_CHAR_MECA(DDL_IMPO=_F(GROUP_NO=("ring_top"),
                                   LIAISON='ENCASTRE'),
                       MODELE=model)


sinkF1 = AFFE_CHAR_MECA(FORCE_NODALE=(_F(GROUP_NO=("bottom_tip"),
                                         FX=0,
                                         FY=0,
                                         FZ=-caseinfo['Frame']['weight_tip'],
                                         ),
                                      _F(GROUP_NO=("ring_bottom"),
                                         FX=0,
                                         FY=0,
                                         FZ=-caseinfo['Frame']['weight_per_metter']*caseinfo['CageShape']['cageCircumference']*9.81/float(meshinfo['NT']),
                                         ),
                                      ),
                        MODELE=model)

listr = DEFI_LIST_REEL(DEBUT=0.0,
                       INTERVALLE=_F(JUSQU_A=tend, PAS=dt))

times = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=listr, PAS_MINI=1e-8),
                       METHODE='AUTO')


for k in range(1, itimes+1):
    INCLUDE(UNITE=91, INFO=0)

# group HDPE pipes
stat2 = CALC_CHAMP(
    RESULTAT=resn,
    GROUP_MA=('topRings',),
    CONTRAINTE=(
        'EFGE_ELGA',
        'EFGE_ELNO',
        'SIPO_ELNO',),
)


# group CABLE
stat3 = CALC_CHAMP(
    RESULTAT=resn,
    GROUP_MA=('twines', 'ropes'),
    CONTRAINTE=(
        'EFGE_ELGA',
        'EFGE_ELNO',
    ),
)

# write results
IMPR_RESU(FORMAT='MED',
          RESU=(_F(CARA_ELEM=elemprop,
                   LIST_INST=listr,
                   NOM_CHAM=('DEPL', 'SIEF_ELGA'),
                   # TOUT_CMP=(DEPL','ACCE','VITE' ),
                   RESULTAT=resn,
                   TOUT_CMP='OUI'),
                _F(GROUP_MA=('topRings',),
                   RESULTAT=stat2,
                   NOM_CHAM=('EFGE_ELGA', 'EFGE_ELNO', 'SIPO_ELNO'
                             ),),
                _F(GROUP_MA=('twines', 'ropes'),
                   RESULTAT=stat3,
                   NOM_CHAM=('EFGE_ELGA', 'EFGE_ELNO',
                             ),),
                ),
          UNITE=80)

# reaction force
stat = CALC_CHAMP(RESULTAT=resn,
                  CONTRAINTE=('SIEF_ELNO',
                              ),
                  FORCE=('REAC_NODA', ),
                  )
# reaction force table
reac1 = POST_RELEVE_T(ACTION=_F(GROUP_NO=('ring_top'),
                                INTITULE='sum reactions',
                                MOMENT=('DRX', 'DRY', 'DRZ'),
                                NOM_CHAM=('REAC_NODA'),
                                OPERATION=('EXTRACTION', ),
                                POINT=(0.0, 0.0, 0.0),
                                RESULTANTE=('DX', 'DY', 'DZ'),
                                RESULTAT=stat))
# write results
IMPR_TABLE(FORMAT_R='1PE12.3',
           TABLE=reac1,
           UNITE=10)

FIN()
