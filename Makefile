all:
	(cd Flux/flux_src; clr; make clean; make)
	(cd Neoclassical/neoclassical_src; clr; make clean; make)
	(cd Phase/phase_src; clr; make clean; make)
	(cd IslandDynamics/islanddynamics_src; clr; make clean; make)
