# UDF scripts for spatially varying porous media model in FLUENT

This collection of UDFs are designed to improve the accuracy of porous media modeling of flow inside coiled aneurysms. It is based on the work by Julia Bhathal who demonstrated that the porosity profile of coiled aneurysms follows a bi-linear behavior. The core or center of the coiled aneurysm has a spatially constant and isotropic porosity, and a region between the core and aneurysm wall has a linear decaying profile (1 at the wall). 


The UDF create bi-linear values of viscous resistance, inertial resistance, and porosity that can be hooked to these parameters inside the PM model of Fluent. 

The functions work by identifying the closest face of the aneurysm dome for each cell, calculating the distance between the face and cell, and applying either a spatially decreasing value or the core value depending on a distance threshold. 

The algorithm can take a few minuted to run upon first initialization as it performs a KD-searching algorithm to find the closest face to each cell in the aneurysm volume.

To do:
* Add porosity function
* test on linux 