all:
	(cd Flux/Source; ../../Scripts/clr; make clean; make)
	(cd Neoclassical/Source; ../../Scripts/clr; make clean; make)
	(cd Phase-1.4/Source; ../../Scripts/clr; make clean; make)
	(cd Phase-2.0/Source; ../../Scripts/clr; make clean; make)
	(cd IslandDynamics/Source; ../../Scripts/clr; make clean; make)
