**ANALYSIS OF YEAST mRNA SEQUENCES**
===================================

*by* Vaskin Maksym

# Introduction

We are going to reuse an existing report, describing a simple analysis protocol, to analyze a new dataset.

## Objectives

* This exercise show many of the strengths of the Unix command-line to process genomic data.
* We will also practice reporting such commands and the results using a text file and the `MarkDown` syntax.

First of all we need to create a project folder, of course, and set it as the current working directory...

```{.sh}
mkdir ex1
cd ex1
```

# Protocol

## The original data

We have to analyze sequence length distribution and GC content for the current mRNA sequences annotated on *S.cerevisiae* genome release Apr. 2011 (SacCer_Apr2011/sacCer3). We have connected to the [UCSC genome browser download web site for yeast](http://hgdownload.soe.ucsc.edu/downloads.html#yeast), and followed the [Full data set link](http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/). We can see that there is a file named `mrna.fa.gz`, we can just copy the corresponding link to the following command:

```{.sh}
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz
```

We obtained a compressed file containing all the mRNA sequences annotated for this genome in `fasta` format. We can check it's content:

```{.sh}
zcat mrna.fa.gz | head -10
#>A18178 1
#caccaataaaaaaacaagcttaacctaattc
#>A21196 1
#cggccagatcta
#>A21197 1
#agcttagatctggccgggg
#>AX557348 1
#gcggatttactcaggggagagcccagataaatggagtctgtgcgtccaca
#gaattcgcacca
#>AX557349 1
```

Let's see how many mRNA sequences do we have. We can uncompress the file and keep going with the flat text file, which can take a lot of disk resources, or we can keep using the compressed file as in the previous command.

```{.sh}
zcat mrna.fa.gz | egrep -c '^>'


zcat mrna.fa.gz | egrep -c '^>' 
#There are 474 mRNAs
```

In case that a program can deal with compressed files, we can take adantage of that feature. There is an `egrep` version that can read such files, if you were guessing is `zegrep` of course (as it happens with `cat` and `zcat`).

```{.sh}
zegrep -c '^>' mrna.fa.gz
#474
```


## Analysis of CDS sequences

To calculate the nucleotide length and the GC content, we can use one of the programs that is provided within the [EMBOSS suite](http://emboss.sourceforge.net/): `infoseq`.

```{.sh}
zcat mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -outfile mrna.lengc.tbl -noheading -only -name -length -pgc
    head -5
#Display basic information about sequences
#A18178         31     29.03  
#A21196         12     58.33  
#A21197         19     63.16  
#AX557348       62     53.23  
#AX557349       62     53.23 
```

In the above command we are just looking to the expected output, the following is doing the job and saving the output into `mrna.lengc.tbl` file.

```{.sh}
zcat mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -outfile mrna.lengc.tbl \
          -noheading -only -name -length -pgc
```

## Visualizing the analysis

By running `R` command, we enter in the `R` shell interpreter, which understands `R` commands

```{.sh}
R
#R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#Copyright (C) 2015 The R Foundation for Statistical Computing
#Platform: x86_64-pc-linux-gnu (64-bit)

#R is free software and comes with ABSOLUTELY NO WARRANTY.
#You are welcome to redistribute it under certain conditions.
#Type 'license()' or 'licence()' for distribution details.

#  Natural language support but running in an English locale

#R is a collaborative project with many contributors.
#Type 'contributors()' for more information and
#'citation()' on how to cite R or R packages in publications.

#Type 'demo()' for some demos, 'help()' for on-line help, or
#'help.start()' for an HTML browser interface to help.
#Type 'q()' to quit R.
```

```{.r}
data <- read.table("mrna.lengc.tbl", header=FALSE)
head(data,4)

#       V1 V2    V3
#1   A18178 31 29.03
#2   A21196 12 58.33
#3   A21197 19 63.16
#4 AX557348 62 53.23
```

Rename the table columns, so they are more meaningful:

```{.r}
colnames(data) <- c("ID","NUClen","GCpct");
head(DATA,4);

#        ID NUClen GCpct
#1   A18178     31 29.03
#2   A21196     12 58.33
#3   A21197     19 63.16
#4 AX557348     62 53.23
```

Let's calculate some stats on the dataset.

```{.r}
summary(data)

#       ID          NUClen           GCpct 
# A18178  :  1   Min.   :  10.0   Min.   :19.00 
# A21196  :  1   1st Qu.:  85.0   1st Qu.:36.83 
# A21197  :  1   Median : 294.5   Median :40.15  
# AB016599:  1   Mean   : 662.0   Mean   :41.90  
# AB109221:  1   3rd Qu.: 897.0   3rd Qu.:46.15  
# AB305041:  1   Max.   :6915.0   Max.   :73.33  
# (Other) :468                                     
```

Now, let's make an histogram:

```{.r}
hist(data$NUClen);
hist(data$GCpct);
```

Or to compare both measures and save into a PNG image...

```{.r}
png(file="plot.png");
plot(data$NUClen ~ data$GCpct);
dev.off();
```

![Showing GC content versus sequence length](plot.png "Showing GC content versus sequence length")


Merging the 3 plots

```{.r}
png(file="plot2.png");
def.par <- par();
# preparing a layout grid where to combine different plots
nf <- layout(matrix(c(2,0,1,3), # matrix contents is plot order
                    2, 2, byrow=TRUE),
             c(3,1), c(1,3), # cell relative sizes
             TRUE);



# computing data distribution
xhist <- hist(data$GCpct, plot=FALSE);
yhist <- hist(data$NUClen, plot=FALSE);
datamax <- max(xhist$counts,
               yhist$counts);

# drawing the main plot
par(mar=c(5,5,1,1)); 
plot(DATA$NUClen ~ DATA$GCpct,
      main="",
      xlab="GC%",
      ylab="Sequence length (bp)",
      col="green");
lines(lowess(data$NUClen ~ data$GCpct),
      col="red", lty="dotted");
mtext(paste("n=",nrow(data),sep=""),
      side=3, line=-1);
# drawing x-axis histogram
par(mar=c(0,5,1,1));
barplot(xhist$counts, ylim=c(0, datamax),
        axes=FALSE, space=0, 
        col="steelblue", main="S.cerevisiae mRNAs")
# drawing y-axis histogram
par(mar=c(5,0,1,1));
barplot(yhist$counts, xlim=c(0, datamax),
        axes=FALSE, horiz=TRUE, space=0, 
        col="steelblue", width=0.5, main="")

par(def.par); # reseting graphical parameters
dev.off();
```


![Showing GC content versus sequence length and marginal distributions](plot2.png "Showing GC content versus sequence length and marginal distributions")


#Analysis per Chromosome

First, RNA information per each chromosome (total of 16 + Mitochondria (M)) was downloaded from (ftp://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/)



```{.sh}
wget --glob --passive-ftp 'ftp://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/database/
chr*_mrna.txt.gz'
```

This data only shows the start and end points of each annotated mRNA, not sequence. We create a new file in which we print the column 11, 15, 10, 17, 18 and 19 of each chromosome file, corresponding, respectively to ID, chromosome, strand, nucleotide start, nucleotide end.

```{.sh}
zcat chr*_mrna.txt.gz | gawk '{ print $11, $15, $10, $17, $18, $19 }' > mrnabychr.tbl
```

Then we create another table from the ID lists of the original mrna.fa.gz, eliminating some FASTA annotations like ">" in front of the ID

```{.sh}
zegrep '^>' mrna.fa.gz | sed 's/^>//; s/ .*$//'> mrna_ids.tbl

#Use this command to obtain pattern from a file

egrep -f mrna_ids.tbl mrnabychr.tbl 

#If we use the word count for the two files created
wc mrnabychr.tbl mrna_ids.tbl
#919  5514 29229 mrnabychr.tbl
#  474   474  3704 mrna_ids.tbl
# 1393  5988 32933 total
```


We can see that there are more annotations in RNAs by chromosomes than IDs in the original file with RNAs of all chromosomes combined. This is because there are copies of same coding genes in different chromosomes and/or in tandem in the same chromosome. To prove that, lets for example count how many times one of the transcripts IDs (K00952) is repeated.

```{.sh}
egrep K00952 mrnabychr.tbl -c
#37
```

Now, merge the information of the table with mRNA lenght and GC percentage with the information of mRNA by chromosome. The problem is that they will not match perfectly because of the problem described above. We overcome that by printing to each ID a numeric extension, starting from 1 and going up +1 for each repetition of the ID.

```{.sh}
gawk 'BEGIN{ OFS="\t"; while(getline<ARGV[1]>0) T[$1]=$2"\t"$3; ARGV[1]=""; } { E[$1]++; print \
$1"."sprintf("%02d",E[$1]), ($1 in T) ? T[$1] : "NA\tNA",$2,$3,$4,$5,$6 }' mrna.lengc.tbl 
mrnabychr.tbl | more 

#We create a new file with this combined info
gawk 'BEGIN{ OFS="\t"; while(getline<ARGV[1]>0) T[$1]=$2"\t"$3; ARGV[1]=""; } { E[$1]++; print \
$1"."sprintf("%02d",E[$1]), ($1 in T) ? T[$1] : "NA\tNA",$2,$3,$4,$5,$6 }' mrna.lengc.tbl 
mrnabychr.tbl > mrna.lengc+chr.tbl 
```

#Creating graphs with R
```{.r}
bychr<-read.table(file="mrna.lengc+chr.tbl",sep="\t",header=FALSE)
summary (bychr)
#V1            V2               V3              V4      V5
#AB016599.01:  1   Min.   :  25.0   Min.   :19.00   chrXII :112   -:464  
#AB016599.02:  1   1st Qu.:  76.0   1st Qu.:38.02   chrVII :108   +:455  
# AB016599.03:  1   Median :  85.0   Median :48.31   chrXIII: 85          
# AB016599.04:  1   Mean   : 385.6   Mean   :46.14   chrIV  : 76          
# AB016599.05:  1   3rd Qu.: 288.5   3rd Qu.:54.55   chrXVI : 73          
# AB016599.06:  1   Max.   :6915.0   Max.   :66.67   chrXV  : 64          
# (Other)    :913                                    (Other):401          
#       V6                V7                V8       
# Min.   :    207   Min.   :    802   Min.   : 1.00  
# 1st Qu.: 194098   1st Qu.: 196400   1st Qu.: 1.00  
# Median : 422936   Median : 423010   Median : 1.00  
# Mean   : 440102   Mean   : 441586   Mean   : 1.59  
# 3rd Qu.: 586635   3rd Qu.: 589664   3rd Qu.: 2.00  
# Max.   :1528679   Max.   :1531783   Max.   :16.00




png(file="plotbychr1.png")
plot(bychr$V2 ~ bychr$V3, col=bychr$V4, main="GC content VS lenght", xlab="GC%", ylab="Lenght")
dev.off();

#Calculating mean mRNA GC% by chromosome and standard deviation
ddd2 <-tapply(bychr$V3, list(bychr$V4), mean)
serror2 <-tapply(bychr$V3, list(bychr$V4), sd)

#Calculating mean mRNA length by chromosome and standard deviation
ddd3 <-tapply(bychr$V2, list(bychr$V4), mean)
serror3 <-tapply(bychr$V2, list(bychr$V4), sd)


```
![Showing GC content VS length for by each chromosome (color-coded by chromosome)](plotbychr1.png "GC% VS mRNA length")

```{.r}
#Making a boxplot for %GC in mRNAs per chromosome
#Standard errors for mRNA length and GC (specially length) are very high, 
#because mRNAs have a very variable length. Outliers were excluded in graphs.

serror3
# chrI      chrII     chrIII      chrIV      chrIX       chrM       chrV 
# 894.44147 1049.16910   30.22878  745.40370  981.15850  193.13791  743.30007 
#     chrVI     chrVII    chrVIII       chrX      chrXI     chrXII    chrXIII 
# 646.99556  616.81148  384.70493  557.80979  681.44858 1006.43484  426.47975 
#    chrXIV      chrXV     chrXVI 
# 528.06859  734.99430  626.73828 

dat <-read.table(file="mrna.lengc+chr.tbl",sep="\t",header=FALSE)
summary (dat)
png(file="boxplot1.png")
dat <- read.table ('mrna.lengc+chr.tbl')
boxplot (dat [dat$V4=="chrI",]$V3, dat [dat$V4=="chrII",]$V3,
		dat [dat$V4=="chrIII",]$V3, dat [dat$V4=="chrIV",]$V3,
		dat [dat$V4=="chrV",]$V3,dat [dat$V4=="chrVI",]$V3,
		dat [dat$V4=="chrVII",]$V3,dat [dat$V4=="chrVIII",]$V3,
		dat [dat$V4=="chrIX",]$V3,dat [dat$V4=="chrX",]$V3,
		dat [dat$V4=="chrXI",]$V3,dat [dat$V4=="chrXII",]$V3,
		dat [dat$V4=="chrXIII",]$V3,dat [dat$V4=="chrXIV",]$V3,
		dat [dat$V4=="chrXV",]$V3,dat [dat$V4=="chrXVI",]$V3,
		dat [dat$V4=="chrM",]$V3, ylim =c(0, 100), outline=FALSE, main="%GC per Chromosome", 
		ylab="%GC", xlab="Chromosomes", names=levels(dat$V4), col="steelblue")
dev.off();




#Making a boxplot for mRNA length per chromosome
png(file="boxplot2.png")
boxplot (dat [dat$V4=="chrI",]$V2, dat [dat$V4=="chrII",]$V2,
		dat [dat$V4=="chrIII",]$V2, dat [dat$V4=="chrIV",]$V2,
		dat [dat$V4=="chrV",]$V2,dat [dat$V4=="chrVI",]$V2,
		dat [dat$V4=="chrVII",]$V2,dat [dat$V4=="chrVIII",]$V2,
		dat [dat$V4=="chrIX",]$V2,dat [dat$V4=="chrX",]$V2,
		dat [dat$V4=="chrXI",]$V2,dat [dat$V4=="chrXII",]$V2,
		dat [dat$V4=="chrXIII",]$V2,dat [dat$V4=="chrXIV",]$V2,
		dat [dat$V4=="chrXV",]$V2,dat [dat$V4=="chrXVI",]$V2,
		dat [dat$V4=="chrM",]$V2, ylim =c(0, 2000), outline=FALSE, main="Length per Chromosome", 
		ylab="Length", xlab="Chromosomes", names=levels(dat$V4), col="steelblue")
dev.off();


```

![Showing GC content by each chromosome](boxplot1.png "GC content per Chromosome")
![Showing GC content by each chromosome](boxplot2.png "mRNA length per Chromosome")



To test if there is any statistically significant differences in GC content and mRNA length, anova was used

```{.r}
a = aov(V3~V4, data=bychr)
summary(a)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#V4           16   9069   566.8   6.128 4.98e-13 ***
#Residuals   902  83429    92.5
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

a2 = aov(V2~V4, data=bychr)
summary(a2)
#             Df    Sum Sq Mean Sq F value  Pr(>F)   
#V4           16  16063162 1003948   2.024 0.00983 **
#Residuals   902 447484158  496102 
```

To see which of the chromosomes are different, Tukey Test was used

```{.r}
tuk<-TukeyHSD(a, conf.level = 0.95)
tuk

#From here, selecting those with p<0.05
#chrVI-chrV        8.73887612   1.0593342 16.41841802 0.0093717
#chrVII-chrVI     -7.21906130 -14.1956248 -0.24249776 0.0337662
#chrXV-chrVI      -9.34923491 -16.8162089 -1.88226097 0.0018850
#chrXVI-chrVI     -7.92681625 -15.2488513 -0.60478119 0.0190334
#chrXV-chrX       -5.95938508 -11.9035601 -0.01521006 0.0486102
#chrM-chrII      -13.31723404 -20.9623095 -5.67215861 0.0000003
#chrM-chrIII     -11.99705882 -22.0083649 -1.98575279 0.0041039
#chrM-chrIV      -11.04263158 -18.0720960 -4.01316711 0.0000087
#chrM-chrIX      -16.07772727 -25.3162535 -6.83920101 0.0000003
#chrV-chrM        10.58629630   3.1446418 18.02795082 0.0001241
#chrVI-chrM       19.32517241  10.7728695 27.87747537 0.0000000
#chrVII-chrM      12.10611111   5.3922972 18.81992500 0.0000001
#chrVIII-chrM     14.97500000   7.3622506 22.58774939 0.0000000
#chrX-chrM        15.93532258   8.6745098 23.19613539 0.0000000
#chrXI-chrM       15.32382353   7.1080271 23.53961991 0.0000000
#chrXII-chrM      11.81107143   5.1247170 18.49742584 0.0000002
#chrXIII-chrM     11.91882353   5.0005038 18.83714321 0.0000004
#chrXIV-chrM      10.28840909   2.5384820 18.03833617 0.0005918
#chrXV-chrM        9.97593750   2.7538493 17.19802570 0.0002479
#chrXVI-chrM      11.39835616   4.3262238 18.47048848 0.0000040


tuk2<-TukeyHSD(a2, conf.level = 0.95)
tuk2

#chrM-chrII      -603.700798 -1163.6039  -43.79770 0.0200885
#chrXIII-chrII   -475.555945  -919.6263  -31.48555 0.0219722
```

We can see that the mRNA GC content varies by chromosome. mRNAs encoded by mitochondrial DNA (M) are the ones with the lowest GC content, reaching statistical significance at 95% confidence level with every other chromosome except chrI. Highest GC content is in mRNAs from chromosome VI, and it reaches statistical significant difference when compared to chromosomes V, VII,XV, XVI and M. There is also  a borderline significant difference between XV and X.

As for mRNA length, they are although averages vary greatly between chromosomes, they are not so statistically significant. Only comparison between ChrII VS Chr XIII and ChrII and mitochondrial mRNA lengthreturns p-values <0.05. Perhaps its due to very large variability of mRNA lengths and/or the fact it doesnt follow normal distribution. Indeed, median of mRNA lengths are less variable.

```{.sh}
#Median mRNA length by chromosome
ddd4 <-tapply(bychr$V2, list(bychr$V4), median)

png(file="mybarplot4.png");
def.par <- par();
par(lwd=2, mar=c(6,4,4,2), mfrow=c(1,1))
mybarplot4<-barplot(ddd4, col="steelblue", main = "Median mRNA length per chromosome", 
xlab="Chromosome number", ylab="Median ofmRNA length", ylim=c(0, 150), cex.main=1.6, cex.sub=1.5, 
cex.names=0.5, space=0.6)
par(def.par)

dev.off();
```

![Showing GC content by each chromosome](mybarplot4.png "mRNA length per Chromosome")





# CODA

## Software

We have used the following versions:

```{.sh}
uname -a
# Linux maksym-VirtualBox 4.10.0-38-generic #42~16.04.1-Ubuntu SMP Tue Oct 10 16:32:20 UTC 2017 
x86_64 x86_64 x86_64 GNU/Linux


wget --version
# GNU Wget 1.17.1 built on linux-gnu.


infoseq -version
# EMBOSS:6.6.0.0

R --version
#R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#Copyright (C) 2015 The R Foundation for Statistical Computing
#Platform: x86_64-pc-linux-gnu (64-bit)
```

##Command used to compile PDF

```{.sh}
pandoc -f markdown_github -t latex      \
       --variable papersize:a4paper     \
       --variable geometry:margin=1.5cm \
       --variable fontsize=11pt         \
       --highlight-style pygments       \
       -o ex2.pdf    \
          ex2.md
```

