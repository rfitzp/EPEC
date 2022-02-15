

EPEC
====

Short Description
-----------------
Suite of programs to simulate island dynamics in presence of RMP

Keywords
--------
tokamak, magnetic island, dynamics, resonant magnetic perturbation (RMP), GPEC code

Long Description
----------------
Suite of programs to simulate multi-harmonic magnetic island chain dynamics 
in presence of resonant magnetic perturbation in time-varying toroidal 
tokamak equilibrium. The simulation scheme is based on asymptotic matching.

Relevant Publications
---------------------
* R. Fitzpatrick, SangKyeun Kim, and Jaehyun Lee.
  *Modeling of q95 windows for the suppression of edge localized modes by resonant magnetic perturbations in the KSTAR tokamak*.
  Phys. Plasmas **28**, 082511 (2021).
* R. Fitzpatrick.
  *Further modeling of q95 windows for the suppression of edge localized modes by resonant magnetic perturbations in the DIII-D tokamak*.
  Phys. Plasmas **28**, 022503 (2021).
* R. Fitzpatrick.
  *Modeling q95 windows for the suppression of edge localized modes by resonant magnetic perturbations in the DIII-D tokamak*.
  Phys. Plasmas **27**, 102511 (2020).
* R. Fitzpatrick, and A.O. Nelson.
  *An improved theory of the response of DIII-D H-mode discharges to static resonant magnetic perturbations and its implications for the suppression of edge localized modes*.
  Phys. Plasmas **27**, 072501 (2020).

Weblinks
--------
* `EPEC Homepage <http://farside.ph.utexas.edu/EPEC-documentation/epec.html>`_
* `Asymptotic Matching Approach to Modeling Tearing Mode Dynamics in Tokamak Plasmas <http://farside.ph.utexas.edu/talks/ifs2021.pdf>`_
* `Predicting Stability Windows in q95 for RMP-Induced ELM Suppression in H-Mode Tokamak Discharges <http://farside.ph.utexas.edu/talks/PoPWebinar.pdf>`_

Programs
--------
/FLUX 
  Program to read gFile(s) and write equilibrium data needed by 
  NEOCLASSICAL and PHASE to fFile(s)

/NEOCLASSICAL 
  Program to read FLUX data from fFile(s), profile data from pFile(s) 
  and cFiles(s), and write neoclassical data needed by PHASE to nFile(s)

/PHASE 
  Program to read FLUX data from fFile(s), NEOCLASSICAL data from 
  nFile(s), GPEC data from uFile(s), mFiles(s), and lFile(s), and 
  perform island dynamics simulation. Final 
  state of plasma saved in sFile.

/RESCALE 
  Program to rescale equilibrium gFile, as well as profile pFile and cFile,
  so as to scan physical quantity as function of time. Allowed scans
  are electron density (n_e), electron temperature (T_e), major radius
  (R_0), toroidal momentum diffusivity (chi_p), total pressure (P),
  ExB frequency (w_E), and edge safety-factor (q_95).
  
/GPEC
  Independent OMFIT module used to calculate ideal-MHD response of tokamak
  plasma equilibirum to resonant magnetic perturbation

Calculation Phases
------------------
* **PHASE 1**: PHASE reads gFile(s) and outputs plasma equilibrium data for Stage 2
* **PHASE 2**: PHASE constructs flux coordinate system, calculates metric quantities, 
  locates resonant surfaces, calculates neoclassical data at resonant surfaces, calculates
  tearing stability matrix, calculates GGJ data, calculates island
  saturation data, and outputs all data to fFile(s). 
* **PHASE 3**: NEOCLASSICAL reads fFile(s), pFile(s), and cFiles(s), calculates neoclassical/profile data at resonant surfaces, 
  and outputs data to nFile(s).
* **PHASE 4**: PHASE reads fFile(s), nFile(s), lfiles(s), uFile(s), mFile(s) for Stages 5 & 6,  and calculates vacuum island widths as function of relative 
  RMP coil phase.
* **PHASE 5**: PHASE performs island dynamics simulation in fixed equilibrium.
* **PHASE 6**: PHASE performs island dynamics simulation in dynamic equilibrium.  
  

Servers
-------
EPEC can be run either directly on the local server (assuming that the EPEC source is installed on the
local server) or via SLURM on a SLURM server, depending on whether
the "SLURM" checkbox on the MAIN screen of the GUI is unchecked or checked, respectively.
There is also a "Batch" checkbox on the MAIN screen: if checked this enables batch execution of
SLURM jobs (when practical and advantageous). The SLURM server can be changed in the MAIN
screen of the GUI. 

Filesystems
-----------
* **Local OMFIT file system.** Usually in /tmp/<username>/OMFIT/ on local machine.
* **Remote OMFIT file system.** Usually in /tmp/<username>/OMFIT/ on SLURM server.
* **Local EPEC-runs database.**
    * In  projectsDir/../EPEC-runs/<Device>.<Shot>.<Time>.<runid>
    * Alternately in projectsDir/EPEC-runs/<Device>.<Shot>.<Times[0]>.<Times[-1]>.<runid> on local machine.
    * Here, <Device>, <Shot>, <Time>, <Times[]>, <runid> are set on the MAIN
      page of the GUI. 
    * Previous runs of EPEC can be loaded from the EPEC-runs database using the
      "Load Previous run ..." button, which will appear on the MAIN page of the GUI if the run is present
      in the database.
    * The present run can be saved to the database by giving "runid" on the
      MAIN page of the GUI a unique value, and then hitting the "Save run"
      button on the MAIN page.

Test Data  
---------
Test data can be loaded into EPEC using the TEST DATA page of the GUI. The following data is available:

* NSTX Shot 127317 at 400 ms, n=1 RMP.
* KSTAR Shot 18594 at 6450 ms, n=2 RMP.
* DIIID Shot 145380 at 3400 ms, n=3 RMP.
* DIIID Shot 145380 between 2000 and 4500 ms, 50 ms interval, n=3 RMP.

Workflows
---------

* **Test run**
  * Hit 'DIIID #145380 3400 ms n=3' button on the TEST DATA page of the GUI.
  * Hit 'Run FLUX + NEOCLASSICAL + PHASE' on the EPEC page of the GUI. 

* **Basic run**
   * Edit the namelists that control FLUX, NEOCLASSICAL, and PHASE:
       * EPEC/INPUTS/epec_namelist - should not need to change.
       * FLUX/INPUTS/flux_namelist - need to choose NTOR, PSIPED.
       * NEOCLASSICAL/INPUTS/neoclassical_namelist - need to choose ExB, NN, LN, YN.
       * PHASE/INPUTS/phase_namelist - need to choose TSTART, TEND, DT, MID, COPT, FREQ.
       * PHASE/INPUTS/rmp_waveform - need to specify RMP waveform. 
   * Load the gFile into FLUX/INPUTS/.
   * Load the pFile into NEOCLASSICAL/INPUTS/. 
   * Load the cFile into NEOCLASSICAL/INPUTS/. 
   * Load the gFile into GPEC.
   * Run GPEC with 1kA (with the correct toroidal mode number) flowing in the desired coil-set. The toroidal mode number 'nn' in the GPEC GUI namelist must match
     NTOR in FLUX/INPUTS/flux_namelist.
   * Hit 'Import l/u/m files from GPEC (attempts to guess coil set)'
     on EPEC page of GUI.
   * Repeat the previous two steps until all coil-sets have been accounted for. (The maximum allowed number of coil-sets is presently 3.)
   * Hit 'Run FLUX + NEOCLASSICAL + PHASE' on the EPEC page of the GUI.
   * The RMP waveform can be modified by changing entries in the PHASE/INPUTS/rmp_waveform namelist. Once the waveform is
     changed, it is only necessary to rerun PHASE (hit the 'Run PHASE' button on the PHASE page of the GUI). 

* **Interpolated run**
   * Edit the namelists that control EPEC, FLUX. NEOCLASSICAL, and PHASE:
       * EPEC/INPUTS/epec_namelist - need to choose TSTART, TEND, DTF, DTP.
       * FLUX/INPUTS/flux_namelist - need to choose NTOR, PSIPED.
       * NEOCLASSICAL/INPUTS/neoclassical_namelist - need to choose ExB, NN, LN, YN.
       * PHASE/INPUTS/phase_namelist - need to choose DT, MID, COPT, FREQ.
       * PHASE/INPUTS/rmp_waveform - need to specify RMP waveform. 
   * Load the gFiles into FLUX/INPUTS/gFiles/. Load the Index into FLUX/INPUTS/gFiles/. The Index should list the g-filenames and the corresponding experimental  times
     in two columns, in order of increasing time.
   * Load the pFiles into NEOCLASSICAL/INPUTS/pFiles/. Load the Index into NEOCLASSICAL/INPUTS/pFiles/. The Index should list
     the p-filenames and the corresponding experimental times in two columns, in order of increasing time.
   * Load the cFiles into NEOCLASSICAL/INPUTS/pFiles/. Load the Index into NEOCLASSICAL/INPUTS/cFiles/. The Index should list
     the c-filenames and the corresponding experimental times in two columns, in order of increasing time.     
   * Load the gFiles into GPEC.
   * Run GPEC with 1kA (with the correct toroidal mode number) flowing in the desired coil-set. The toroidal mode number 'nn' in GPEC GUI must match
     NTOR in FLUX/INPUTS/flux_namelist.
   * Hit 'Import l/u/m files from GPEC (attempts to guess coil set)'
     on the EPEC page of the GUI.
   * Repeat previous two steps until all coil-sets have been accounted for. (The maximum allowed number of coil-sets is presently 3.)
   * Hit 'Generate fFiles and nFiles' on the EPEC page of the GUI.
   * Hit 'Run EPEC' on the EPEC page of the GUI.
   * The RMP waveform can be modified by changing entries in the PHASE/INPUTS/rmp_waveform namelist. Once the waveform is
     changed, you can hit 'Run EPEC' again.

* **Rescaled run**
   * Edit the namelists that control EPEC, RESCALE, FLUX. NEOCLASSICAL, and PHASE:
       * EPEC/INPUTS/epec_namelist  - need to choose TSTART, TEND, DTF, DTP and RESCALE parameters. 
       * FLUX/INPUTS/flux_namelist - need to choose NTOR, PSIPED.
       * NEOCLASSICAL/INPUTS/neoclassical_namelist - need to choose ExB, NN, LN, YN.
       * PHASE/INPUTS/phase_namelist - need to choose DT, MID, COPT, FREQ.
       * PHASE/INPUTS/rmp_waveform - need to specify RMP waveform. 
   * Load the gFile into FLUX/INPUTS/.
   * Load the pFile into NEOCLASSICAL/INPUTS/. 
   * Load the cFile into NEOCLASSICAL/INPUTS/.
   * Hit one of the buttons on the SCAN page of the GUI. This will rescale the gFile/pFile/cFile (e.g., by changing q_95).
     The rescaled files will be loaded into FLUX/INPUTS/gFiles, NEOCLASSSICAL/INPUTS/pFiles, and NEOCLASSICAL/INPUTS/cFiles,
     together with the appropriate Index files. 
   * Load the gFiles into GPEC.
   * Run GPEC with 1kA (with the correct toroidal mode number) flowing in the desired coil-set. The toroidal mode number 'nn' in GPEC GUI must match
     NTOR in FLUX/INPUTS/flux_namelist.
   * Hit 'Import l/u/m files from GPEC (attempts to guess coil set)'
     on the EPEC page of the GUI.
   * Repeat the previous two steps until all coil-sets have been accounted for. (The maximum allowed number of coil-sets is presently 3.)
   * Hit 'Generate fFiles and nFiles' on the EPEC page of the GUI.
   * Hit 'Run EPEC' on the EPEC page of the GUI.
   * The RMP waveform can be modified by changing entries in the PHASE/INPUTS/rmp_waveform namelist. Once the waveform is
     changed, you can hit 'Run EPEC' again. 
    
    
