#updated on Jan 8,2024 for stratified sampling
#Sep 21, 2023 
#Updated on Nov 29 2023
#APPLES needs same size of sequences.  So for the CO1 gene I added hyphen to the begining and end
#EPA-NG querry should have the same size as reference. 
library(foreach)
#install.packages("phytools")
library(phytools)
library(maps)
library(phangorn)
library(phylotools)
library(dplyr)
getwd()
#Read fasta files of placed co1 on the bb tree
R_co1_20_1=read.fasta(file = "CO1_20_1S.fasta")#This fasta has 80%co1, the bb tree is 20% #named as bb percentage for co1
R_co1_20_2=read.fasta(file = "CO1_20_2S.fasta")
R_co1_20_3=read.fasta(file = "CO1_20_3S.fasta")
R_co1_20_4=read.fasta(file = "CO1_20_4S.fasta")
R_co1_20_5=read.fasta(file = "CO1_20_5S.fasta")
R_co1_20_6=read.fasta(file = "CO1_20_6S.fasta")
R_co1_20_7=read.fasta(file = "CO1_20_7S.fasta")
R_co1_20_8=read.fasta(file = "CO1_20_8S.fasta")
R_co1_20_9=read.fasta(file = "CO1_20_9S.fasta")
R_co1_20_10=read.fasta(file = "CO1_20_10S.fasta")
R_co1_40_1=read.fasta(file = "CO1_40_1S.fasta")
R_co1_40_2=read.fasta(file = "CO1_40_2S.fasta")
R_co1_40_3=read.fasta(file = "CO1_40_3S.fasta")
R_co1_40_4=read.fasta(file = "CO1_40_4S.fasta")
R_co1_40_5=read.fasta(file = "CO1_40_5S.fasta")
R_co1_40_6=read.fasta(file = "CO1_40_6S.fasta")
R_co1_40_7=read.fasta(file = "CO1_40_7S.fasta")
R_co1_40_8=read.fasta(file = "CO1_40_8S.fasta")
R_co1_40_9=read.fasta(file = "CO1_40_9S.fasta")
R_co1_40_10=read.fasta(file = "CO1_40_10S.fasta")
R_co1_60_1=read.fasta(file = "CO1_60_1S.fasta")
R_co1_60_2=read.fasta(file = "CO1_60_2S.fasta")
R_co1_60_3=read.fasta(file = "CO1_60_3S.fasta")
R_co1_60_4=read.fasta(file = "CO1_60_4S.fasta")
R_co1_60_5=read.fasta(file = "CO1_60_5S.fasta")
R_co1_60_6=read.fasta(file = "CO1_60_6S.fasta")
R_co1_60_7=read.fasta(file = "CO1_60_7S.fasta")
R_co1_60_8=read.fasta(file = "CO1_60_8S.fasta")
R_co1_60_9=read.fasta(file = "CO1_60_9S.fasta")
R_co1_60_10=read.fasta(file = "CO1_60_10S.fasta")
R_co1_80_1=read.fasta(file = "CO1_80_1S.fasta")
R_co1_80_2=read.fasta(file = "CO1_80_2S.fasta")
R_co1_80_3=read.fasta(file = "CO1_80_3S.fasta")
R_co1_80_4=read.fasta(file = "CO1_80_4S.fasta")
R_co1_80_5=read.fasta(file = "CO1_80_5S.fasta")
R_co1_80_6=read.fasta(file = "CO1_80_6S.fasta")
R_co1_80_7=read.fasta(file = "CO1_80_7S.fasta")
R_co1_80_8=read.fasta(file = "CO1_80_8S.fasta")
R_co1_80_9=read.fasta(file = "CO1_80_9S.fasta")
R_co1_80_10=read.fasta(file = "CO1_80_10S.fasta")
R_co1_99_1=read.fasta(file = "CO1_99_1S.fasta")
R_co1_99_2=read.fasta(file = "CO1_99_2S.fasta")
R_co1_99_3=read.fasta(file = "CO1_99_3S.fasta")
R_co1_99_4=read.fasta(file = "CO1_99_4S.fasta")
R_co1_99_5=read.fasta(file = "CO1_99_5S.fasta")
R_co1_99_6=read.fasta(file = "CO1_99_6S.fasta")
R_co1_99_7=read.fasta(file = "CO1_99_7S.fasta")
R_co1_99_8=read.fasta(file = "CO1_99_8S.fasta")
R_co1_99_9=read.fasta(file = "CO1_99_9S.fasta")
R_co1_99_10=read.fasta(file = "CO1_99_10S.fasta")

#check the number of characters in swequence column 
nchar(R_co1_99_1$seq.text)
#check the class
class(R_co1_80_1)
dim(R_co1_99_1)
#20% tree has 3283 co1, 40% has 2469, 60% has 1654 co1, 80% has 813, 99% has 27 co1
#Total 24809 (pre 2709, co1 682, post 21418)
# pre_dashes <- data.frame(rep(strrep("-",2709),3616))
# post_dashes <- data.frame(rep(strrep("-",21418),3616))
pre_dashes <- data.frame(rep(strrep("-",2709),27))
post_dashes <- data.frame(rep(strrep("-",21418),27))

#post_dashes <- data.frame(rep(strrep("-",24127),3616))
colnames(pre_dashes) <- "preSeq"
colnames(post_dashes) <- "postSeq"
21418+2709+682



# class(pre_dashes$rep. <- strrep......681...904.)
# aa <- coalesce(pre_dashes$rep.strrep......681...904.,R_co1_20_1$seq.text,post_dashes$rep.strrep......21417...904.)
# aa <- data.frame(aa)
# nchar(aa$aa)    
df_R_co1_99_1 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_1$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_2 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_2$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_3 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_3$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_4 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_4$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_5 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_5$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_6 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_6$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_7 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_7$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_8 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_8$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_9 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_9$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_10 <- data.frame(paste(pre_dashes$preSeq,R_co1_99_10$seq.text,post_dashes$postSeq, sep = ""))

pre_dashes <- data.frame(rep(strrep("-",2709),3283))
post_dashes <- data.frame(rep(strrep("-",21418),3283))
colnames(pre_dashes) <- "preSeq"
colnames(post_dashes) <- "postSeq"

df_R_co1_99_1 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_1$seq.text,post_dashes$postSeq, sep = ""))
#df_R_co1_80_1 <- data.frame(paste(R_co1_80_1$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_2 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_2$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_3 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_3$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_4 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_4$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_5 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_5$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_6 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_6$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_7 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_7$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_8 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_8$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_9 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_9$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_10 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_10$seq.text,post_dashes$postSeq, sep = ""))

pre_dashes <- data.frame(rep(strrep("-",2709),828))
post_dashes <- data.frame(rep(strrep("-",21418),828))
colnames(pre_dashes) <- "preSeq"
colnames(post_dashes) <- "postSeq"

df_R_co1_99_1 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_1$seq.text,post_dashes$postSeq, sep = ""))
#df_R_co1_80_1 <- data.frame(paste(R_co1_80_1$seq.text,pos2_dashes$postSeq, sep = ""))
df_R_co1_99_2 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_2$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_3 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_3$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_4 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_4$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_5 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_5$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_6 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_6$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_7 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_7$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_8 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_8$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_9 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_9$seq.text,post_dashes$postSeq, sep = ""))
df_R_co1_99_10 <- data.frame(paste(pre_dashes$preSeq,R_co1_20_10$seq.text,post_dashes$postSeq, sep = ""))



colnames(df_R_co1_99_1) <- "seq"
colnames(df_R_co1_99_2) <- "seq"
colnames(df_R_co1_99_3) <- "seq"
colnames(df_R_co1_99_4) <- "seq"
colnames(df_R_co1_99_5) <- "seq"
colnames(df_R_co1_99_6) <- "seq"
colnames(df_R_co1_99_7) <- "seq"
colnames(df_R_co1_99_8) <- "seq"
colnames(df_R_co1_99_9) <- "seq"
colnames(df_R_co1_99_10) <- "seq"

# colnames(df_R_co1_80_1) <- "seq"
# colnames(df_R_co1_80_2) <- "seq"
# colnames(df_R_co1_80_3) <- "seq"
# colnames(df_R_co1_80_4) <- "seq"
# colnames(df_R_co1_80_5) <- "seq"
# colnames(df_R_co1_80_6) <- "seq"
# colnames(df_R_co1_80_7) <- "seq"
# colnames(df_R_co1_80_8) <- "seq"
# colnames(df_R_co1_80_9) <- "seq"
# colnames(df_R_co1_80_10) <- "seq"


df_R_co1_99_1_new <- data.frame(cbind(R_co1_20_1$seq.name,df_R_co1_99_1$seq))
df_R_co1_99_2_new <- data.frame(cbind(R_co1_20_2$seq.name,df_R_co1_99_2$seq))
df_R_co1_99_3_new <- data.frame(cbind(R_co1_20_3$seq.name,df_R_co1_99_3$seq))
df_R_co1_99_4_new <- data.frame(cbind(R_co1_20_4$seq.name,df_R_co1_99_4$seq))
df_R_co1_99_5_new <- data.frame(cbind(R_co1_20_5$seq.name,df_R_co1_99_5$seq))
df_R_co1_99_6_new <- data.frame(cbind(R_co1_20_6$seq.name,df_R_co1_99_6$seq))
df_R_co1_99_7_new <- data.frame(cbind(R_co1_20_7$seq.name,df_R_co1_99_7$seq))
df_R_co1_99_8_new <- data.frame(cbind(R_co1_20_8$seq.name,df_R_co1_99_8$seq))
df_R_co1_99_9_new <- data.frame(cbind(R_co1_20_9$seq.name,df_R_co1_99_9$seq))
df_R_co1_99_10_new <- data.frame(cbind(R_co1_20_10$seq.name,df_R_co1_99_10$seq))

class(df_R_co1_80_1_new)
nchar(df_R_co1_99_1_new$X2)
#Generate fasta files for new co1 dataset
library(seqRFLP)
dataframe2fas(df_R_co1_99_1_new, file = "df_S_co1_20_1_new.fasta")
dataframe2fas(df_R_co1_99_2_new, file = "df_S_co1_20_2_new.fasta")
dataframe2fas(df_R_co1_99_3_new, file = "df_S_co1_20_3_new.fasta")
dataframe2fas(df_R_co1_99_4_new, file = "df_S_co1_20_4_new.fasta")
dataframe2fas(df_R_co1_99_5_new, file = "df_S_co1_20_5_new.fasta")
dataframe2fas(df_R_co1_99_6_new, file = "df_S_co1_20_6_new.fasta")
dataframe2fas(df_R_co1_99_7_new, file = "df_S_co1_20_7_new.fasta")
dataframe2fas(df_R_co1_99_8_new, file = "df_S_co1_20_8_new.fasta")
dataframe2fas(df_R_co1_99_9_new, file = "df_S_co1_20_9_new.fasta")
dataframe2fas(df_R_co1_99_10_new, file = "df_S_co1_20_10_new.fasta")
getwd()
