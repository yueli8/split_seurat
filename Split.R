library(Seurat) #加载seurat软件
library(stringr)
library(tidyverse)
a=Read10X("D:/GSE178318/All") #读取3个文件
a <- CreateSeuratObject(a) #创建Seurat对象
head(a@meta.data) #查看读取的数据前6行，结果如下展示
head(a@meta.data, 8) #查看前8行的命令
#正式拆分，使用tidyverse函数
library(tidyverse) #加载tidyverse
rownames(a@meta.data) -> a@meta.data$rowLeo #反向赋值，增加了1列
head(a@meta.data) #展示如下
str_split(a$rowLeo,"_")  #将-作为拆分的字符，拆的结果如下展示
#展示前6行
head(str_split(a$rowLeo,"_",simplify=T)) #展示前6行
head(str_split(a$rowLeo,"_",simplify=T) [,2])  #展示第2列的前6行
str_split(a$rowLeo,"_",simplify=T) [,2] -> a@meta.data$Sample #反向赋值，实现分组
str_split(a$rowLeo,"_",simplify=T) [,3] -> a@meta.data$Tissue #反向赋值
head(a@meta.data)
tail(a@meta.data)
#按照样本（Sample）将数据提取出来
Leo <- SplitObject(a,split.by='Sample')
Leo
#按组织（Tissue）提取
Xuan <- SplitObject(a,split.by='Tissue') 
Xuan
#将seurat对象提取出来，并保存
head(Xuan$CRC@meta.data)
head(Xuan$CRC@meta.data$Sample,100000) 
table(Xuan$CRC@meta.data$Sample)
head(Xuan$CRC@meta.data$Sample,17001)
Xuan2 <- SplitObject(Xuan$CRC,split.by='Sample') 
Xuan2
Xuan3 <- SplitObject(Xuan$LM,split.by='Sample') 
Xuan3
Xuan4 <- SplitObject(Xuan$PBMC,split.by='Sample') 
Xuan4
CRC07=Xuan2$COL07
CRC12=Xuan2$COL12
CRC15=Xuan2$COL15
CRC16=Xuan2$COL16
CRC17=Xuan2$COL17
CRC18=Xuan2$COL18
LM07=Xuan3$COL07
LM12=Xuan3$COL12
LM15=Xuan3$COL15
LM16=Xuan3$COL16
LM17=Xuan3$COL17
LM18=Xuan3$COL18
PBMC12=Xuan4$COL12
PBMC17=Xuan4$COL17
PBMC18=Xuan4$COL18
saveRDS(CRC07, file = "CRC07.Rds")
saveRDS(CRC12, file = "CRC12.Rds")
saveRDS(CRC15, file = "CRC15.Rds")
saveRDS(CRC16, file = "CRC16.Rds")
saveRDS(CRC17, file = "CRC17.Rds")
saveRDS(CRC18, file = "CRC18.Rds")

saveRDS(LM07, file = "LM07.Rds")
saveRDS(LM12, file = "LM12.Rds")
saveRDS(LM15, file = "LM15.Rds")
saveRDS(LM16, file = "LM16.Rds")
saveRDS(LM17, file = "LM17.Rds")
saveRDS(LM18, file = "LM18.Rds")

saveRDS(PBMC12, file = "PBMC12.Rds")
saveRDS(PBMC17, file = "PBMC17.Rds")
saveRDS(PBMC18, file = "PBMC18.Rds")