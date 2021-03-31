# RevPublish

Steps to process the data and perform analyses:

ROOT: the project root directory.

(A) To create an aggregated databases for the multi-unit activities (MUAs) and
(single-unit activities (SUs), run:

	(a) ROOT/ANALYZE_RAW/MAIN_AGGREGATE_MUA_DATA.M, and b)
	(ROOT/ANALYZE_RAW/MAIN_AGGREGATE_SU_DATA.M 
	
	This will create 2 data files, one for the MUA and another for SU, in the
	ROOT/_DATA DIRECTORY.
	
(B) Analyze: (1) Use ANALYZE_RECONSTRUCT/CREATE_BF_TABLE.M to create a file
with the best-frequencies (BFs). This procedure creates a table of BFs that is
saved in the _DATA directory. This new file ends with _BF. The software will
use this file to sort neurons according to their "quality". Creates a
DATA_MUA_XXX_BF.MAT file.

	(2) Run ANALYZE_RECONSTRUCT/MAIN_LOOPOVER_UNITS.M to perform the
	(reconstruction on various number of units. Each such case is save seperately
	(at _DATA/_RECONSTRUC. Creates the RECONSTRUCT_(MUA)_XXX.MAT or
	(RECONSTRUCT_(MUA)_XXX.MAT files.

	(3) Run ANALYZE_RAW/MAIN_RAW_MUA_STATS.M to get statistics on RAW
	(measurements, before averaging. This will creates the CCSTATS_(MUA)_XXX.MAT files.
