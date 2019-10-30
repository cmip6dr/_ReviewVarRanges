# _ReviewVarRanges
Review Variable Ranges

The aim is to review ranges (extremes and percentiles), comparing across models to detect outliers of unexpected behaviour. The software should produce 2 to 4 plots for each variable.

For CMIP5, scanning was done over the historical experiment. Ideally it would be done over a broader range of experiments. It may make sense to group experiments, e.g. provide one set of values for:
(a) historical and control experiments (historical, piControl);
(b) storyline scenarios (1pctCO2, abrupt4xCO2, sspxxx);
(c) sensitivity studies;

The percentiles are estimated, based on percentiles in each time slot or on a sample of time slots. This allows analysis to be computed on a file by file basis, and then aggregated across files.


For each time slot:
(a) 2x5 extreme values, 2xcount of extreme value occurence, 2x10 at 10-4, 2x10 at 10-3, and 98 percentile values (150 values).
(b) Aggregating: extreme of extremes, aggregate extreme value counts. Median percentiles?


T-Digest: https://arxiv.org/pdf/1902.04023.pdf appears to be a useful concept -- more systematic than approach used for CMIP5.

Workflow: (1) identify files and time slices, (2) decide on strategy (all data vs sampled time slices), (3) step through.

(1) file names, # time slots per file, data shape.
(2) Percentiles and extremes -- for each variable & experiment group.


