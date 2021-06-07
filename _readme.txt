Steps to process the data and perform analyses:

ROOT: the project root directory.

(A) To create an aggregated databases for the multi-unit activities (MUAs) and
(single-unit activities (SUs), run:

	(a) ROOT/ANALYZE_RAW/MAIN_AGGREGATE_MUA_DATA.M, and b)
	(ROOT/ANALYZE_RAW/MAIN_AGGREGATE_SU_DATA.M 
	
	This will create 2 data files, one for the MUA and another for SU, in the
	ROOT/_DATA DIRECTORY.
	
(B) Analyze: 
	(1) ANALYZE_RECONSTRUCT/CREATE_BF_TABLE.M: use it to create a file
	with the best-frequencies (BFs). This procedure creates a table of BFs that is
	saved in the _DATA directory. This new file ends with _BF. The software will
	use this file to sort neurons according to their "quality". Creates a
	DATA_MUA_XXX_BF.MAT file.

	(2) ANALYZE_RECONSTRUCT/MAIN_LOOPOVER_UNITS.M: use it to perform the
	(reconstruction on various number of units. Each such case is save seperately
	(at _DATA/_RECONSTRUC. Creates the RECONSTRUCT_(MUA)_XXX.MAT or
	(RECONSTRUCT_(MUA)_XXX.MAT files.

	(3) ANALYZE_RAW/MAIN_RAW_MUA_STATS.M: use it to get statistics on RAW
	(measurements, before averaging. This will creates the CCSTATS_(MUA)_XXX.MAT files.

	(4) ANALYZE_FEATURES/MAIN_CREATE_CC_DATABASE.M: compares spectrogram reconstructions (estimations) with other DRR conditions and saves the database into ANALYZED_CC_XX mat file.





(C) Plots:
	(1) _FIGS4PAPES/RESULTS_CC_ERRORBAR_ANALYSIS.M: plots 2 bar plots with erro bars, one for SU and another for MUA. In each figure the bars compare spectrogram reconstructions of Sdry-to-Sest (blue bars) and Sdrr-to-Sest (red bars). 
	
	(2) _FIGS4PAPES/RESULTS_CC_ANALYSIS.M: 
		(a) Two plots of SU & MUA reconstructions. Each plot shows CC as a function of DRRs for various number of units used for the reconstructions.
		(b) Boxplot of SU & MUA side by side on the same figure.
		(c) Two-ways ANOVA plot (data-type vs DRRs).




(*?*) Analysis & plot results for the paper:  -	Run
ANALYZE_RECONSTRUCT/ANALYZE_RECONSTRUCTION.M: this will analyze -	the
spectrogram 	 reconstructions between the DRY and other DRR conditions and
will plot all relevant figures.

-	Run ANALYZE_RAW/ANALYZE_STATS.M to analyze for: * RMD: response modulation
depth * MG : modulation gain * CCr: correlation coefficient between DRY and
other DRR conditions

-	Run ANALYZE_FEATURES\MAIN_CREATE_CC_DATABASE.M to create a database file
@_DATA/Analysis. This  database contains CC between various DRR along the time
domain.
	
-	Use ANALYZE_FEATURES\ANALYZE_CC.M to analyze the CCs and to plot results.



