# UiS-Aqua
--- A hydrodynamic module for aquaculture structures.

**The main aim for the development of this module is to provide a controllable and repeatable numerical solution for reliable analyses of aquaculture structures.** 


> Although considerable numerical programs have been proposed with a large number of publications, most of the proposed numerical program are either  commercial or in-house. In order to promote the availability  and accessibility of numerical simulations for aquaculture structures, a hydrodynamic modules is developed here for conducting time-domain structural analysis of aquaculture structures with [code aster](https://www.code-aster.org/V2/spip.php?rubrique2). 



### Work environments

In order to use this module, you need to install the following packages or programs first:

* [python 3.6 (or higher)](https://www.python.org/)

* [Numpy 1.21 (or higher)](https://numpy.org/)

* [Code_aster 14.6 (or higher stable version)](https://www.code-aster.org/spip.php?article272)

  Please check if `as_run` located at the following default installation path:

  > /opt/aster146/bin/as_run

  If not, you need to change Line 21 in [`run.sh`](./Example/run.sh) and Line 2 in [`ASTERRUN.export`](./Example/asterinput/ASTERRUN.export) according to your environments.

* [salome_meca 2019.0.3 (or higher)](https://www.code-aster.org/spip.php?article303)

  Please check if `salome` located at the following default installation path:

  Default installation path:

  > /opt/salome_meca/appli_V2019.0.3_universal/salome

  If not, you need to change Line 7 in [`run.sh`](./Example/run.sh) according to your environments.



### Run your first case

The [Example](./Example/README.MD) folder provide a simple setup for a squared fish cage under pure current conditions. 

1. Open a terminal and change the directory to the `Example`

2. Type the command:

   ``` shell
   sh run.sh
   ```

3. Wait until the job finish. If it shows `exit_code=0` at the end, it means the simulation finish without any error.

```shell
 --------------------------------------------------------------------------------
 Copying results

copying .../fort.80...                                                  [  OK  ]
copying .../fort.10...                                                  [  OK  ]
copying .../fort.6...                                                   [  OK  ]

<A>_ALARM          Code_Aster run ended


 
 ---------------------------------------------------------------------------------
                                            cpu     system    cpu+sys    elapsed
 ---------------------------------------------------------------------------------
   Preparation of environment              0.00       0.00       0.00       0.00
   Copying datas                           0.05       0.02       0.07       0.08
   Code_Aster run                        823.59      22.92     846.51     812.87
   Copying results                         0.01       0.03       0.04       0.04
 ---------------------------------------------------------------------------------
   Total                                 823.87      23.04     846.91     813.28
 ---------------------------------------------------------------------------------

as_run 2020.0

------------------------------------------------------------
--- DIAGNOSTIC JOB : <A>_ALARM
------------------------------------------------------------


EXIT_CODE=0
Change the variable back to the default value....>>>>>>

```

4. Want to clean the generated file and run a different cases? use `sh clean.sh ` to clean these files and run again. 

   

### Document

* [Case setup](./doc/demon.md)
* [Mesh generator](./doc/mesh.md)
* [Hydrodynamic module for line-type elements](doc/1d_hydro.md)
* [Hydrodynamic module for plane-type elements](doc/2d_hydro.md)
* [Wake effects](doc/wakeEffect.md)
* [Waves](doc/waves.md)

### Publications 

  1. Cheng, H., Li, L., Aarsæther, K. G., & Ong, M. C. (2020). Typical hydrodynamic models for aquaculture nets: A comparative study under pure current conditions. Aquacultural Engineering, 90, 102070. https://doi.org/10.1016/j.aquaeng.2020.102070

  2. Cheng, H., Ong, M. C., Li, L., & Chen, H. (2022). Development of a coupling algorithm for fluid-structure interaction analysis of flexible nettings in fluid. Ocean Engineering, 243, 110208. https://doi.org/10.1016/j.oceaneng.2021.110208

  3. Mjåtveit, M. A., Cheng, H., Ong, M. C., & Lee, J. (Under review). Numerical study of two typical gravity-based fish cages with different dimensions under pure current conditions. Aquacultural Engineering.

  4. Cheng, H., Li, L., & Ong, M. C. (Under review). Structural responses of gravity type fish cages under pure current conditions. Aquaculture.

