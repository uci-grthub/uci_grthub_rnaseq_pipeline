#!/usr/bin/env Rscript 

library("iSEEde")
library("airway")
library("DESeq2")
library("iSEE")

# Example data ----

data("airway")
airway$dex <- relevel(airway$dex, "untrt")

dds <- DESeqDataSet(airway, ~ 0 + dex + cell)

dds <- DESeq(dds)
res_deseq2 <- results(dds, contrast = list("dextrt", "dexuntrt"))
head(res_deseq2)
#> log2 fold change (MLE): dextrt vs dexuntrt 
#> Wald test p-value: dextrt vs dexuntrt 
#> DataFrame with 6 rows and 6 columns
#>                   baseMean log2FoldChange     lfcSE      stat      pvalue       padj
#>                  <numeric>      <numeric> <numeric> <numeric>   <numeric>  <numeric>
#> ENSG00000000003 708.602170     -0.3812539  0.100654 -3.787752 0.000152016 0.00128292
#> ENSG00000000005   0.000000             NA        NA        NA          NA         NA
#> ENSG00000000419 520.297901      0.2068127  0.112219  1.842944 0.065337213 0.19646961
#> ENSG00000000457 237.163037      0.0379205  0.143445  0.264356 0.791505314 0.91141884
#> ENSG00000000460  57.932633     -0.0881679  0.287142 -0.307054 0.758802543 0.89500551
#> ENSG00000000938   0.318098     -1.3782416  3.499906 -0.393794 0.693733216         NA

# iSEE / iSEEde ---

airway <- embedContrastResults(res_deseq2, airway, name = "dex: trt vs untrt")
contrastResults(airway)
#> DataFrame with 63677 rows and 1 column
#>                   dex: trt vs untrt
#>                 <iSEEDESeq2Results>
#> ENSG00000000003 <iSEEDESeq2Results>
#> ENSG00000000005 <iSEEDESeq2Results>
#> ENSG00000000419 <iSEEDESeq2Results>
#> ENSG00000000457 <iSEEDESeq2Results>
#> ENSG00000000460 <iSEEDESeq2Results>
#> ...                             ...
#> ENSG00000273489 <iSEEDESeq2Results>
#> ENSG00000273490 <iSEEDESeq2Results>
#> ENSG00000273491 <iSEEDESeq2Results>
#> ENSG00000273492 <iSEEDESeq2Results>
#> ENSG00000273493 <iSEEDESeq2Results>

saveRDS(airway, "airway.rds")
