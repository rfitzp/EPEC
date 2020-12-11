all:
	(cd Rescale/Source; ../../Scripts/clr; make clean; make)
	(cd Flux/Source; ../../Scripts/clr; make clean; make)
	(cd Neoclassical/Source; ../../Scripts/clr; make clean; make)
#	(cd Phase-1.x/Source; ../../Scripts/clr; make clean; make)
	(cd Phase-2.x/Source; ../../Scripts/clr; make clean; make)
	(cd IslandDynamics/Source; ../../Scripts/clr; make clean; make)
