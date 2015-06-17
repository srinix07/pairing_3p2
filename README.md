# pairing_3p2
3p2 pairing khodel.

## Prerequisite
* gsl libraries in standard location.
* make
* g++(gcc) 
* OpenMP (Optional)

## Usage
* Run make inside khodel_3p2 directory.
* Then **./a.out --help**

## NN-Potential
* *Entem-Machleidt* (500Mev) potential is given in _pot_ directory.
* The format of potential file is, 
* k (fm^-1) k' (fm^-1) v11 v13 v31 v33 (Mevfm^3)

## Example Usage
* **./a.out -p pot/test_vsrg_N3LO_EM_500_NN_3P2-3F2_Mev_fm3.dat > gap.dat**
* To get gap in Mev,
*  **./a.out -p pot/test_vsrg_N3LO_EM_500_NN_3P2-3F2_Mev_fm3.dat -M > gap.dat**
* For relativistic calculation use *-r* flag.
* See *./a.out --help* for other options.

## OpenMP.
* By default Openmp is used, to disable OpenMp,
* remove -fopenmp flags from *CFLAGS* and *LIB* in the makefile.
