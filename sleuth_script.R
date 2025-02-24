library('sleuth') # load sleuth library
stab = read.table("samples.txt",header=TRUE) # read table
so = sleuth_prep(stab) # prepare sleuth object using table

so = sleuth_fit(so, ~condition, 'full') # fit full model with condition as explanatory variable
so = sleuth_fit(so, ~1, 'reduced') # fit a reduced model
so = sleuth_lrt(so, 'reduced', 'full') # likelihood ratio test to compare reduced and full model

library(dplyr) # load dplyr library

sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) # get the results from test above

# filter significant results with an adjusted p-value (q-value) threshold of 0.05
# sort the results in ascending order of raw p-values
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 
head(sleuth_significant, n = 10) # view only first ten

# write the table to a text file with only target_id, pval, and qval and top 10 results
write.table(sleuth_significant, file="sleuth_results.txt",quote = FALSE,row.names = FALSE)
head(dplyr::select(sleuth_significant, target_id, pval, qval), n=10)
