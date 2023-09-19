# CAD-O

> ...
## Usage
### Files upload online and import to QGIS
Choose the geographical zone (tile, canton) for the area of study. Step 1 : Download the files, Step 2 : Import them to QGIS

Swissbuildings3D : 
1. Geodatabase format from https://www.swisstopo.admin.ch/en/geodata/landscape/buildings3d3.html#download
2. Import *Floor/Roof/Wall* layers

MO cadaster
1. GeoPackage format from https://geodienste.ch/services/av
2. Import *lcsf* layer, filter by "Genre"="batiment"
### QGIS
Create a new GeoPackage file. Add layers with the "export features" option of QGIS (export without "fid" field), with names :
- *zone_cad* : MO features of buildings connected to DHN
- *zone_tout* : MO features of all buildings in the area of study
- *centrale* : point feature of thermal heating station coordinates
- *floor* : swissbuildings3D features of floors of all buildings in the area of study
- *wall* : swissbuildings3D features of walls of all buildings in the area of study
- *roof* : swissbuildings3D features of roofs of all buildings in the area of study

### Code custom modifications
In addition to the newly created GeoPackage, a climatic and a horizon file must be provided and imported on the same directory as where the software code was cloned.

Before launching the "main_code.py" from command line, it must be custom modified for each simulation, either from a text files reading app or a python editor app such as *Spyder* or *Visual Code Studio*.
Each modification is signaled with *#TODO*.
- gpkg_filepath = r"---.gpkg" : File path of the GeoPackage containing all necessary layers
- create_geometry_3D = True/False (default = False) : Activates the simulation with thermal envelope from Swissbuildings3D geometries (much longer simulation)
- calculate_volume_3D = True/False (default = False) : Activates the volume calculation from Swissbuildings3D geometries
- citysim_filepath = r"---/CitySim.exe" : File path of the CitySim solver
- directory_path = r"---" : Name of the new directory to be created by the simulation
- climate_file = r"---.cli" : File path of the climate file
- horizon_file = r"---.hor" : File path of the horizon file

### Python libraries
The required libraries to import are listed in the *requirements.txt* file
