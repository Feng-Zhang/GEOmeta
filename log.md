### to do
- [] 测试数据GSE38900太大，不方便测试。选择小一点的数据。
- [x] Undefined global functions or variables:    GPLlist
- [x] 使用fread读取数据。
- [x] feat_20211206_name: `write.table()`默认输出行名，但是行名那一列没有列名。而要用fread的话，需要手动加一列ID。
- [x] fix-20211202-convert-name：修复saveGSE报Error in tapply错误。
- [x] metaGene中pheGene的列名被限定死了，能不能改？如果不能改在Rd文档中进行说明。pheGene的对象要预先限定，要不然data.table的对象传进来会报错。
- [x] saveGSE使用write.table时有些基因名会变成1-Dec字样。这是由于Table(getGEO(GPL,getGPL=TRUE))获取GPL信息时就已经变成这样，没办法修正。
- [x] perf_20211215_uniqueFeature: uniqueFeature只接受data.table数据，而不接受tbl数据，需要把uniqueFeature改成S3对象，接受data.frame, data.table, tbl对象
- [x] feat-20220714: 保持代码风格一致；

