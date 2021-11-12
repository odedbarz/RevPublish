Steps to process the data and perform analyses:

ROOT: the project root directory.

(A) To create an aggregated databases for the multi-unit activities (MUAs) and
(single-unit activities (SUs), run:

	(a) ROOT/ANALYZE_RAW/MAIN_AGGREGATE_MUA_DATA.M, and
	(b) ROOT/ANALYZE_RAW/MAIN_AGGREGATE_SU_DATA.M

	This will create 2 data files, one for the MUA and another for SU, in the
	ROOT/_DATA DIRECTORY.


(B) Analyze:
	(1) ANALYZE_RECONSTRUCT/CREATE_BF_TABLE.M:
		* Use it to create a file with the best-frequencies (BFs). 
		* This procedure creates a table of BFs that is	saved in the _DATA directory. 
		* This new file ends with _BF. 
		* The software will use this file to sort neurons according to their "quality". Creates a	DATA_MUA_XXX_BF.MAT file.

	(2) ANALYZE_RECONSTRUCT/MAIN_LOOPOVER_UNITS.M:
		* Use it to perform the	reconstruction on various number of units. 
		*Each such case is save seperately (at _DATA/_RECONSTRUC. Creates the RECONSTRUCT_(MUA)_XXX.MAT or(RECONSTRUCT_(MUA)_XXX.MAT files.

	(3) ANALYZE_RAW/MAIN_RAW_MUA_STATS.M:
		* Use it to get statistics on RAW measurements, before averaging. 
		* This will creates the CCSTATS_(MUA)_XXX.MAT files.

	(4) ANALYZE_RECONSTRUCT/MAIN_CREATE_CC_DATABASE.M:
		* Compares spectrogram reconstructions (estimations) with other DRR conditions and saves the database into ANALYZED_CC_XX mat file.
		* Creates the CCt (temporal correlation coefficients).

	(5) ANALYZE_RECONSTRUCT/MAIN_LOOPOVER_STATS.m: used to check effect of dependence between reconstructions.
		The analysis is done in _FIGS4PAPER\ADDITIONALS.m. This file uses data from,
		_DATA\STATS\CCUNITS_(09-JUN-2021)_FOR_VALID_STATS_TEST.m.
		
	MULTI-RUNs to compare between BEF sorting and random sorting
		(6) ANALYZE_RECONSTRUCT/MAIN_LOOPOVER_UNITS_MULTIRUNS.M: 
			* Repeats analysis for a given number of units; 
			* The units are selected at random in each run;
			* The results are saved at _DATA/RECONSTRUCT/.
		
		(7) ANALYZE_RECONSTRUCT/MAIN_CREATE_CC_DATABASE_MULTIRUNS.M:
			* Analyze all the runs into one cell and saves it at the _DATA\ANALYSIS folder.

(C) Plots:
	(1) FIGS4PAPER/RESULTS - SINGLE UNIT - CC & RMD MG/MAIN_ANALYZE_CCRESPONSE.M:
		Plots the CC between each response (SU or MU) and the "best envelope". The best envelope is taken from CREATE_BF_TABLE.M (see above). It's the frequency with the highest Pearson correlation (CC) between one of the spectrogram frequency bands and the response. 

	(2) _Figs4Paper/RESULTS - POPULATION RESPONSE/RESULTS_CC_ANALYSIS.M:
		(a) Two plots of SU & MUA reconstructions. Each plot shows CC as a function of DRRs for various number
			of units used for the reconstructions.
			- SUvsMUA_CC(Sdry-vs-Sest)_vs_drr_units(10-100).png
			- CC(Sdry-vs-Sest & Sdrr-vs-Sest)_for_SU_&_MUA_units(100).png
			- SUvsMUA_CC(Sdry-vs-Sest)_vs_drr_units(100).png
		(b) Boxplot of SU & MUA side by side on the same figure.
		(c) Wilcoxon signed rank tests

	(3) _FIGS4PAPER/RESULTS - CCT/MAIN_PLOT_CCT.M:
		(a) Voice vs. Unvoiced histograms.
		(b) 


	(?) _Figs4Paper/RESULTS_CC_ERRORBAR_ANALYSIS.M: plots 2 bar plots with error bars, one for SU and another for MUA. In each figure the bars compare spectrogram reconstructions of Sdry-to-Sest (blue bars) and Sdrr-to-Sest (red bars).










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
