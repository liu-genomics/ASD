library(data.table)
spidex_zscore_file <- "//media/yuwen/Elements/spidex_database/spidex_public_noncommercial_v1_0.tab.dpsi_zscore"
spidex_zscore <- fread(spidex_zscore_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# note that there is Inf and -Inf in the dpsi_zscore
cutoff <- quantile(spidex_zscore$dpsi_zscore, seq(0,1,0.1))
