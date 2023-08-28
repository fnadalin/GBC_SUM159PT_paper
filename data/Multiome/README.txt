cistopic_input_peaks_1C_ARC.bed : all input peaks for cisTopic that match some 1E peaks (adjacent peaks, d <= 0; original coordinates on 1C)

cistopic_repr_regions_non_redundant.bed : all regions marked as "reproducible" across topic pairs 1C-1E (topics with high correlation with sequencing depth are excluded; union of matching peaks coordinates in 1C and 1E according to idr)

MODULE_1.bed : reproducible regions for module 1, 1C-1E. column 5 is min(IDR,1000), columns 6-7 are cisTopic probabilities in 1C and 1E (union of matching peaks coordinates in 1C and 1E according to idr)


