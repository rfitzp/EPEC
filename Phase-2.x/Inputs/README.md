# Input files for program PHASE

## Contents

- Phase.in: Fortran_90 namelist file
- Waveform.in: Fortran_90 namelist file specifying RMP coil current waveform
- fFile: fFile (usually symbolic link)
- /fFiles: Directory containing fFiles used for interpolated equilibria (usually symbolic link)
- nFile: nFile (usually symbolic link)
- /nFiles: Directory containing nFiles used for interpolated equilibria (usually symbolic link)
- uFile: uFile (usually symbolic link)
- /uFiles: Directory containing uFiles used for interpolated equilibria (usually symbolic link)
- mFile: mFile (usually symbolic link)
- /mFiles: Directory containing mFiles used for interpolated equilibria (usually symbolic link)
- lFile: lFile (usually symbolic link)
- /lFiles: Directory containing lFiles used for interpolated equilibria (usually symbolic link)

*/fFiles must contain nFiles as well as Index which lists fFile names and experimental times in two columns*

*/nFiles must contain nFiles as well as Index which lists nFile names and experimental times in two columns*

*/uFiles must contain nFiles as well as Index which lists uFile names and experimental times in two columns*

*/mFiles must contain mFiles as well as Index which lists mFile names and experimental times in two columns*

*/lFiles must contain lFiles as well as Index which lists lFile names and experimental times in two columns*
