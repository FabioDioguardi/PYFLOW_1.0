Name of the program: PYFLOW
Title of the Manuscript: PYFLOW: a computer code for the calculation of the impact parameters of Dilute Pyroclastic Density Currents (DPDC) basing on field data
Authors: Fabio Dioguardi(1), Pierfrancesco Dellino(1)
(1): Dipartimento di Scienze della Terra e Geoambientali, Università degli Studi di Bari “Aldo Moro”, Bari, Italy

PYFLOW is a computer code designed for quantifying the hazard related to Dilute Pyroclastic Density Currents (DPDC). 
DPDCs are multiphase flows that form during explosive volcanic eruptions. They are the major source of hazard related to volcanic eruptions, 
as they exert a significant stress over buildings and transport significant amounts of volcanic ash, which is hot and unbreathable. 
The program calculates the DPDC’s impact parameters (e.g. dynamic pressure and particle volumetric concentration) and is founded on the 
turbulent boundary layer theory adapted to a multiphase framework. Fluid-dynamic variables are searched with a probabilistic approach, 
meaning that for each variable the average, maximum and minimum solutions are calculated. From these values, PYFLOW creates probability 
functions that allow to calculate the parameter at a given percentile. The code is written in Fortran 90 and can be compiled and installed 
on Windows, Mac OS X, Linux operating systems (OS). A User’s manual is provided, explaining the details of the theoretical background, 
the setup and running procedure and the input data. The model inputs are DPDC deposits data, e.g. particle grainsize, layer thickness, 
particles shape factor and density. PYFLOW reads input data from a specifically designed input file or from the user’s direct typing by command lines. 
Guidelines for writing input data are also contained in the package. PYFLOW guides the user at each step of execution, asking for additional data and inputs. 
The program is a tool for DPDC hazard assessment and, as an example, an application to the DPDC deposits of the Agnano-Monte Spina eruption (4.1 ky BP) 
at Campi Flegrei (Italy) is presented.     