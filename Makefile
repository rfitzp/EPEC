all:
	(cd Flux/Source; ../../Scripts/clr; make clean; make)
	(cd Neoclassical/Source; ../../Scripts/clr; make clean; make)
	(cd Phase/Source; ../../Scripts/clr; make clean; make)
	(cd IslandDynamics/Source; ../../Scripts/clr; make clean; make)
