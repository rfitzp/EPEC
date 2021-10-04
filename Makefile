all:
	(cd Flux/Source; ../../Scripts/clr; make clean; make)
	(cd Neoclassical/Source; ../../Scripts/clr; make clean; make)
	(cd Phase/Source; ../../Scripts/clr; make clean; make)
	(cd Rescale/Source; ../../Scripts/clr; make clean; make)
	(cd Sphinx; make html)
