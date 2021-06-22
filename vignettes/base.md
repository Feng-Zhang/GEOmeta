# overview

对GEOquery包进行封装，方便下载多个芯片数据，然后用封装的meta包进行meta分析。

## 对芯片其中一个基因进行meta分析

待定

## 对芯片的所有基因进行meta分析

### 1. 下载所有芯片的表型和基因型

首先，需要知道GSE的id号。如果要做meta分析，还需要知道GSE里哪些个体，即GSM号，和分组信息。我们以内置数据`phe_test`为例：

    data("phe_test")
    gseIds=unique(phe_test$GSE)

如果这步不想做meta分析，只是想下载数据，可以直接提供GSE的id号：

    gseIds=c("GSE32138","GSE24132","GSE38900")

然后进行下载：

    for(i in gseIds){
      print(i)
      GSEtoExpr(i,destdir = "../tmp")
    }
    #> [1] "GSE32138"
    #> [1] "GSE24132"
    #> [1] "GSE38900"

GEOmeta会自动下载所有芯片数据并存储在当前工作目录的tmp文件夹下。

### 2. 合并所有芯片成一个大矩阵

我们需要把下载的芯片的表达矩阵进行合并

    exprs = combineExprs(destdir = "../tmp",GSEs = gseIds)
    dim(exprs)
    #> [1] 46811   269

### 3. 进行meta分析

由于第2步把芯片下所有个体都进行了合并，有时候我们并不是对所有个体都感兴趣，所以我们需要对整合后的矩阵进行过滤。

    data("phe_test")
    inds = intersect(phe_test$GSM,colnames(exprs))
    phe_test = phe_test[phe_test$GSM %in% inds,]
    exprs = unique(exprs[,inds])
    dim(exprs)
    #> [1] 43269    32

接着再对每个基因进行meta分析

    allMeta = NULL
    genes = row.names(exprs)
    for(i in 1:20){
      #i=1
      pheGene = cbind(phe_test,t(exprs[1,]))
      metaDat = metaGene(pheGene,genes[i])
      res = metaRes(metaDat)
      allMeta = rbind(allMeta,res)
    }
    head(allMeta)
    #>   n.e  mean.e     sd.e n.c   mean.c     sd.c  studlab         TE      seTE     lower     upper
    #> 1   6 1.58206 2.435823   6 2.550582 5.499047 GSE38900 -0.2102189 0.5797196 -1.346448 0.9260106
    #> 2   4     NaN       NA   4      NaN       NA GSE32138        NaN        NA       NaN        NA
    #> 3   6     NaN       NA   6      NaN       NA GSE24132        NaN        NA       NaN        NA
    #> 4   6     NaN       NA   6      NaN       NA GSE38900        NaN        NA       NaN        NA
    #> 5   4     NaN       NA   4      NaN       NA GSE32138        NaN        NA       NaN        NA
    #> 6   6     NaN       NA   6      NaN       NA GSE24132        NaN        NA       NaN        NA
    #>         zval  statistic      pval  w.fixed w.random   fixed_TE fixed_lower fixed_upper    fixed_z
    #> 1 -0.3626217 -0.3626217 0.7168875 2.975528 2.975528 -0.2102189   -1.346448   0.9260106 -0.3626217
    #> 2        NaN        NaN       NaN 0.000000 0.000000         NA          NA          NA         NA
    #> 3        NaN        NaN       NaN 0.000000 0.000000         NA          NA          NA         NA
    #> 4        NaN        NaN       NaN 0.000000 0.000000         NA          NA          NA         NA
    #> 5        NaN        NaN       NaN 0.000000 0.000000         NA          NA          NA         NA
    #> 6        NaN        NaN       NaN 0.000000 0.000000         NA          NA          NA         NA
    #>   fixed_pvalue  random_TE random_lower random_upper   random_z random_pvalue I2 tao2 heter_pvalue
    #> 1    0.7168875 -0.2102189    -1.346448    0.9260106 -0.3626217     0.7168875 NA   NA           NA
    #> 2           NA         NA           NA           NA         NA            NA NA   NA           NA
    #> 3           NA         NA           NA           NA         NA            NA NA   NA           NA
    #> 4           NA         NA           NA           NA         NA            NA NA   NA           NA
    #> 5           NA         NA           NA           NA         NA            NA NA   NA           NA
    #> 6           NA         NA           NA           NA         NA            NA NA   NA           NA

### to do

-   测试数据GSE38900太大，不方便测试。选择小一点的数据。  
-   使用fread读取数据。`write.table()`默认输出行名，但是行名那一列没有列名。如果后面要用fread的话，需要手动加一列symbolID，再输出，或者用`write.csv`输出csv文件。
-   combineExprs中的GSEs参数
-   Undefined global functions or variables: GPLlist
