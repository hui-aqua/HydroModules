echo "Prepare....>>>>>>" 
cd asterinput
ln -s ../../src module
cd ..

# generate mesh
/opt/salome_meca/appli_V2019.0.3_universal/salome -t  ./asterinput/module/mesh/SSM.py

# prepare files for code aster
echo "Prepare file >>>>>>"
mkdir asteroutput
mkdir pythonOutput
mv *.med asterinput/inputMesh.med
VAR=$(shuf -i 555-1000 -n 1)
sed -i "s|PATH_OF_THE_CURRENT_FOLDER|$(pwd)|g" asterinput/asterInput1.py
sed -i "s|PATH_OF_THE_CURRENT_FOLDER|$(pwd)|g" asterinput/ASTERRUN.export
sed -i "s|RANDOMNUMBER|$VAR|g" asterinput/ASTERRUN.export

# run the simulation using code aster
echo "Run....>>>>>>"
/opt/aster146/bin/as_run asterinput/ASTERRUN.export 
# 2> $(pwd)/error.log

echo "Change the variable back to the default value....>>>>>>"
sed -i "s|$VAR|RANDOMNUMBER|g" asterinput/ASTERRUN.export
sed -i "s|$(pwd)|PATH_OF_THE_CURRENT_FOLDER|g" asterinput/ASTERRUN.export
sed -i "s|$(pwd)|PATH_OF_THE_CURRENT_FOLDER|g" asterinput/asterInput1.py
