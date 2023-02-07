library(DESeq2)
library(stringr)

# file_list
txt_files = list.files("input_txt/", pattern=".txt")

# set each parameters (replicates, FDR, Figure)
N = 3
FDR = 0.1
fig = c(400, 380)

# create directly
dir.create("result")

for (i in txt_files) {
    # set input
    input = paste0("input_txt/", i)
    output_txt = paste0("result/", str_sub(i, end = -5), "_DESeq2.txt")
    output_png = paste0("result/", str_sub(i, end = -5), "_DESeq2.png")

    print(paste0("Input  txt: ", input))
    print(paste0("Output txt: ", output_txt))
    print(paste0("Output png: ", output_png))

    # read txt
    count = read.table(input, header = TRUE, row.names = 1, sep = "\t", quote="")
    # print(head(count))

    # set group
    group = data.frame(condition = as.factor(c(rep(1, N), rep(2, N))))

    # deseq2
    dds = DESeqDataSetFromMatrix(countData = round(count), colData = group, design = ~ condition) # カウントデータは正の整数である必要があるのでround()で丸める
    # dds$condition = relevel(dds$condition, ref = A) # negativee control 
    dds = DESeq(dds)
    tmp = results(dds)

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

    # FDR閾値(q.value < FDR)を満たす遺伝子数を表示
    print(paste0("FDR < 0.1 : ", sum(q.value < FDR)))    
    # sum(p.adjust(p.value, method="BH") < FDR)

    # table
    result = cbind(rownames(count), count, p.value, q.value, ranking, log2FC, lfcSE, baseMean, stat)

    # output
    write.table(result, file = output_txt, sep = "\t", append = FALSE, quote = FALSE, row.names = FALSE)
    png(output_png, pointsize=14, width=500, height=400)
    plotMA(dds)

}