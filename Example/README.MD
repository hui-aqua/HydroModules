## A squared fish cage under pure current conditions

The following figure is the results of the example case, showing the deformations of a fish cage due to current forces. 

![](/home/hui/GitHub/HydroModules/doc/figures/demon_deformation.png)

### The file structures:

* asterinput
  * module:  A soft link point to the ```./src``` folder, in which the source code are stored. 
  * asterinput1.py: The calculation script, including materials properties, time schemes, output files, etc.
  * asterinput2.py: The calculation script, including numerical schemes, hydrodynamic forces assignment, etc.  
  * ASTERRUN.export --> submit job for code aster, including system environments, version of code aster, etc. 
* run.sh: Run the simulation
* clean.sh: Clean the generated files
* setting.json: Parameter for the case setting. 
* asteroutput:  Generated folder for the out put of med file. 
* pythonOutput: Generated folder for the out put of time-series results.

