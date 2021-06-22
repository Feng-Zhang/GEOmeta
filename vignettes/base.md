# overview

对GEOquery包进行封装，方便下载多个芯片数据，然后用封装的meta包进行meta分析。
\#\# 对芯片其中一个基因进行meta分析 待定 \#\#
对芯片的所有基因进行meta分析 \#\#\# 1. 下载所有芯片的表型和基因型
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

    exprs = combineExprs(destdir = "../tmp")

### 3. 进行meta分析

由于第2步把芯片下所有个体都进行了合并，有时候我们并不是对所有个体都感兴趣，所以我们需要对整合后的矩阵进行过滤。

    data("phe_test")
    inds = intersect(phe_test$GSM,colnames(exprs))
    phe_test = phe_test[phe_test$GSM %in% inds,]
    exprs = unique(exprs[,inds])

接着再对每个基因进行meta分析

    allMeta = NULL
    genes = row.names(exprs)
    for(i in 1:20){
      #i=1
      print(i)
      pheGene = cbind(phe_test,t(exprs[1,]))
      metaDat = metaGene(pheGene,genes[i])
      res = metaRes(metaDat)
      allMeta = rbind(allMeta,res)
    }
    #> [1] 1
    #> [1] 2
    #> [1] 3
    #> [1] 4
    #> [1] 5
    #> [1] 6
    #> [1] 7
    #> [1] 8
    #> [1] 9
    #> [1] 10
    #> [1] 11
    #> [1] 12
    #> [1] 13
    #> [1] 14
    #> [1] 15
    #> [1] 16
    #> [1] 17
    #> [1] 18
    #> [1] 19
    #> [1] 20

### to do

-   测试数据GSE38900太大，不方便测试。选择小一点的数据。  
-   使用fread读取数据。`write.table()`默认输出行名，但是行名那一列没有列名。如果后面要用fread的话，需要手动加一列symbolID，再输出，或者用`write.csv`输出csv文件。
-   combineExprs中的GSEs参数
-   Undefined global functions or variables: GPLlist
