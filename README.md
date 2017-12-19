# gsca_data_process
give a process for data used in gsca project
## cnv
1. percantage calculation use threshold data.
cnv are classified into 4 types, homo amp, homo deletion, hete amp, hete deletion. percentage is calculated.

2. cnv correlation with gene expression use cnv raw data.
use lm(exp ~ cnv) to get the relationship between gene expression and cnv.
fdr is filtered with threshold 0.05. 

## snv
1. percentage
per=mut count/sample count

2. maf
filter raw maf data, only some rows are remained.
cancer type infomation is added as clinical part of maf, and maf part to make cancer type selection available.
save them as maf object.

3. survival
survival difference is calculated between mut and non mut genes.
only genes with > 5% mut frequency will be calculated.(samles > 50)
only genes with > 3 mut frequency will be calculated.(samles < 50)

coxp and logrank p are all get.
data with no fdr filtered.

## methylation
1. methylation difference
diff= tumor(methy) - normal(methy)
p value is calculate by t.test.
fdr = -log10(fdr)
fdr = ifelse(fdr > 50, 50, fdr)
filtered by filter(fdr < 0.05)
direction = estimate > 0 ~ "Down", # normal high
        estimate < 0 ~ "Up" # tumor high
when diff methy is get, cancers with less than 10 normal-tumor paired sample will be droped, and only 14 cancers are remained with diff data.


2. methy correlate with gene expression
use lm(exp ~ cnv) to get the relationship between gene expression and cnv.
fdr is filtered with threshold 0.05.

3. methylation survival difference
strtified use middle methylation level
Hyper_worse = ifelse(estimate > 0,  "High","Low")
## ...
