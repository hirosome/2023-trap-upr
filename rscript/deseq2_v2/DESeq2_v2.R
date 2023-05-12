library(DESeq2)
library(stringr)

# sample ID
sample_ID = c("Total_DMSOvsTPG", "P-TRAP_DMSOvsTPG", "L10a-TRAP_DMSOvsTPG", "P0vsL10a_DMSO")
sample_ID

# set each parameters (replicates, FDR)
N = 3
FDR = 0.1

# set group
group = data.frame(condition = as.factor(c(rep(1, N), rep(2, N))))

# create directly
dir.create("result")

for (input in sample_ID) {
  
  # 全リードの解析
  # read total readCount file
  input_readCount = paste0("input_txt/", input, "_All.txt")
  output_txt = paste0("result/", str_sub(input), "_All_DESeq2.txt")
  print(paste0("Input  txt: ", input_readCount))
  print(paste0("Output txt: ", output_txt))
  
  readCount = read.table(input_readCount, header=TRUE, row.names=1, sep ="\t", quote="")
  head(readCount)
  
  # DESeq2 for total readCount
  dds_readCount = DESeqDataSetFromMatrix(countData=readCount, colData=group, design=~condition)
  dds_readCount = estimateSizeFactors(dds_readCount)
  dds_readCount = estimateDispersions(dds_readCount)
  dds_readCount = nbinomWaldTest(dds_readCount)
  tmp = results(dds_readCount)
  head(tmp)
  
  # 結果を格納
  p.value = tmp$pvalue
  p.value[is.na(p.value)] = 1 # NAを1に置換
  q.value = tmp$padj
  q.value[is.na(q.value)] = 1
  ranking = rank(p.value)
  log2FC = tmp$log2FoldChange
  lfcSE = tmp$lfcSE
  baseMean = tmp$baseMean
  stat = tmp$stat # Wald検定の結果
  
  print(paste0("FDR < 0.1 : ", sum(q.value < FDR))) 
  result = cbind(rownames(readCount), readCount, p.value, q.value, ranking, log2FC, lfcSE, baseMean, stat)
  write.table(result, file=output_txt, sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)
  
  # T>C変換リードの解析
  # read tcReadCount file
  input_tcReadCount = paste0("input_txt/", input, "_TC.txt")
  output_txt = paste0("result/", str_sub(input), "_TC_DESeq2.txt")
  print(paste0("Input  txt: ", input_tcReadCount))
  print(paste0("Output txt: ", output_txt))
  
  tcReadCount = read.table(input_tcReadCount, header=TRUE, row.names=1, sep ="\t", quote="")
  head(tcReadCount) 
  
  # DESeq2 for TcReadCount
  dds_tcReadCount = DESeqDataSetFromMatrix(countData=tcReadCount, colData=group, design=~condition)
  # dds_tcReadCount = estimateSizeFactors(dds_tcReadCount)     # TcReadCountのsize factorは使わない
  sizeFactors(dds_tcReadCount) = sizeFactors(dds_readCount)    # 代わりにreadCountのsize factorを使う
  dds_tcReadCount = estimateDispersions(dds_tcReadCount)
  dds_tcReadCount = nbinomWaldTest(dds_tcReadCount)
  tmp = results(dds_tcReadCount)
  head(tmp)
  
  p.value = tmp$pvalue
  p.value[is.na(p.value)] = 1 # NAを1に置換
  q.value = tmp$padj
  q.value[is.na(q.value)] = 1 # NAを1に置換
  ranking = rank(p.value)
  log2FC = tmp$log2FoldChange
  lfcSE = tmp$lfcSE
  baseMean = tmp$baseMean
  stat = tmp$stat # Wald検定の結果
  
  print(paste0("FDR < 0.1 : ", sum(q.value < FDR))) 
  result = cbind(rownames(tcReadCount), tcReadCount, p.value, q.value, ranking, log2FC, lfcSE, baseMean, stat)
  write.table(result, file = output_txt, sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE)
  
  
  # Non-T>Cリードの解析
  # read nonTcReadCount file
  input_nonTcReadCount = paste0("input_txt/", input, "_nonTC.txt")
  output_txt = paste0("result/", str_sub(input), "_nonTC_DESeq2.txt")
  print(paste0("Input  txt: ", input_nonTcReadCount))
  print(paste0("Output txt: ", output_txt))
  
  nonTcReadCount = read.table(input_nonTcReadCount, header=TRUE, row.names=1, sep ="\t", quote="")
  head(nonTcReadCount) 
  
  # DESeq2 for nonTcReadCount
  dds_nonTcReadCount = DESeqDataSetFromMatrix(countData=nonTcReadCount, colData=group, design=~condition)
  # dds_tcReadCount = estimateSizeFactors(dds_tcReadCount)     # TcReadCountのsize factorは使わない
  sizeFactors(dds_nonTcReadCount) = sizeFactors(dds_readCount)    # 代わりにreadCountのsize factorを使う
  dds_nonTcReadCount = estimateDispersions(dds_nonTcReadCount)
  dds_nonTcReadCount = nbinomWaldTest(dds_nonTcReadCount)
  tmp = results(dds_nonTcReadCount)
  head(tmp)
  
  p.value = tmp$pvalue
  p.value[is.na(p.value)] = 1 # NAを1に置換
  q.value = tmp$padj
  q.value[is.na(q.value)] = 1 # NAを1に置換
  ranking = rank(p.value)
  log2FC = tmp$log2FoldChange
  lfcSE = tmp$lfcSE
  baseMean = tmp$baseMean
  stat = tmp$stat # Wald検定の結果
  
  print(paste0("FDR < 0.1 : ", sum(q.value < FDR))) 
  result = cbind(rownames(nonTcReadCount), tcReadCount, p.value, q.value, ranking, log2FC, lfcSE, baseMean, stat)
  write.table(result, file = output_txt, sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE)
  
}

