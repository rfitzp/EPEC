all:
	(cd Rescale/Source; ../../Scripts/clr; make clean; make)
	(cd Flux/Source; ../../Scripts/clr; make clean; make)
	(cd Neoclassical/Source; ../../Scripts/clr; make clean; make)
	(cd Phase/Source; ../../Scripts/clr; make clean; make)
	(cd IslandDynamics/Source; ../../Scripts/clr; make clean; make)
	(cd fFileGenerate/Source; ../../Scripts/clr; make clean; make)
	(cd WindowFind/Source; ../../Scripts/clr; make clean; make)

