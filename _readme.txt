Steps to process the data and perform analyses:

ROOT: the project root directory.

(1) To create an aggregated databases for the multi-unit activities (MUAs) and 
    single-unit activities (SUs), run:

	(a) ROOT/ANALYZE_RAW/MAIN_AGGREGATE_MUA_DATA.M 
	(b) and ROOT/ANALYZE_RAW/MAIN_AGGREGATE_MUA_DATA.M 
	
	This will create 2 data files, one for the MUA and another for SU, in the 
    ROOT/_DATA DIRECTORY.
	
(2) Use ANALYZE_RECONSTRUCT/BEST_ENVELOPE_FREQUENCY.M to create a file with the 
    best-frequencies (BFs). This procedure creates a table of BFs that is saved in the 
    _DATA directory. This new file ends with _BF. The software will use this file to sort
    neurons according to their "quality".

(3) Run the ANALYZE_RECONSTRUCT/MAIN_LOOPOVER_UNITS.M to perform the reconstruction on various
	number of units. Each such case is save seperately at _DATA/_RECONSTRUC.  
