library('sleuth')
stab = read.table("samples.txt",header=TRUE)
so = sleuth_prep(stab)

so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

library(dplyr)

sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval)
head(sleuth_significant, n = 10)

write.table(sleuth_significant, file="sleuth_results.txt",quote = FALSE,row.names = FALSE)
head(dplyr::select(sleuth_significant, target_id, pval, qval), n=10)
