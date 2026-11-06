Perfmorming scans of selection sweep signatures along the genome for different populations using the tool SweepFinder2 (SF2), and comparing them statistically with PicMin.

**Run SF2**:
1. subsetVCF_batch1234.sh
2. filter_parseVCF_sweepfinder_poparray.sh
3. sf2_input_poparray.sh
4. globalsfs_poparray.sh
5. sf2_pol_chromarray_A1.sh [runs for one population at a time, faster]
6. sf2_pol_poparray.sh [runs for many populations at a time, slower]

**Visualisation SF2**:
1. manhplot_methodcomparison.R [...]
2. manhplot_A1_allscaf.R
3. manhplot_allp_alls.R

**Run PicMin**:
1. picmin_prep.R
2. picmin_onescaf.R
3. picmin_chromarray.sh [runs picmin_onescaf.R for all chromosomes]

**Visualisation PicMin**:
1. picmin_plot_allscafs.R

**Combine plots**:
1. combine_pngplots.R
