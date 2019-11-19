#ディレクトリの準備とツールのインストール
　
$ cd ~　
$ mkdir rrbs
$ cd rrbs
$ mkdir rawdata tools ref fastqc trim map

$ brew install bowtie2
$ brew install samtools

$ cd ~/rrbs/tools
$ wget https://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.22.1.tar.gz
$ tar xvzf bismark_v0.22.1.tar.gz  
$ mv Bismark_v0.22.1/* /usr/local/bin/
$ bismark -v

$ R  #ターミナル上でRを起動する。以下Rコマンド。

> install.packages("BiocManager")

#初期状態ならばインストールに失敗するはず。
#XQuartzの配布サイト(https://www.xquartz.org)からXQuartzをダウンロード＆インストールする。

> install.packages("BiocManager")
> BiocManager::install("methylKit")
> BiocManager::install("genomation")

> packageVersion("BiocManager")
> packageVersion("methylKit")
> packageVersion("genomation")


> q() #Rを終了するコマンド

#Rコマンドここまで。　

#解析データダウンロード(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66121)
#使用するデータ
SRR1812639　血液サンプル　細胞種: CD19+ B-cells  慢性リンパ性白血病患者
SRR1812671　血液サンプル　細胞種: CD19+ B-cells  健常者
 
$ cd ~/rrbs/rawdata 
$ fastq-dump --gzip SRR1812671
$ fastq-dump --gzip SRR1812639

$ ls -lh　　

$ gunzip  SRR1812639.fastq.gz 
$ gunzip  SRR1812671.fastq.gz 

$ ls -lh


$ cd ~/rrbs/ref/ 

#viを起動、010_download-ucsc.shという新しいファイルを作成する。 
$ vi 010_download-ucsc.sh     
       

#「i」を入力しmインサートモードへ切り替え

~
~
-- INSERT --

#インサートモードの状態で下記のように入力。

#!/bin/bash
set -euo pipefail
u1="ftp://hgdownload.soe.ucsc.edu"
u2="goldenPath/hg38/bigZips/analysisSet"
u3="hg38.analysisSet.chroms.tar.gz"
curl -O ${u1}/${u2}/${u3}
tar zxvf ${u3}

#この後に「:wq」と入力し、010_download-ucsc.shファイルに入力した内容を保存。

$ less 010_download-ucsc.sh

#中身を確認したら「q」で元の画面に戻る。

$ chmod a+x 010_download-ucsc.sh
$ ./010_download-ucsc.sh
 

$ ls

$ for i in $(seq 1 22)
> do
> mv hg38.analysisSet.chroms/chr${i}.fa  ~/rrbs/ref
> done

$ ls


#遺伝子アノテーションファイルダウンロード
#UCSC Genome Browserから取得する。まずUCSC Genome Browsr(https://genome.ucsc.edu/index.html)を開く。
#設定画面では下記のように設定する。設定を確認したら「get output」をクリックし、次の画面に移動しよう。
  
　assembly: Dec.2013(GRCh38/hg38)
　group: Genes and Gene Predictions
　track: NCBI RefSeq
　table: UCSC RefSeq (refGene)
　region: genome
　output format: custom track
　output file: refGene.bed
　file type returned: plain text


$ mv ~/Downloads/refGene.bed ~/rrbs/ref
$ ls 


#情報解析
#FastQCで確認

$ cd ~/rrbs
$ fastqc --nogroup -o fastqc/ -t 4 rawdata/*.fastq

$ open fastqc/SRR1812639_fastqc.html 
$ open fastqc/SRR1812671_fastqc.html

#トリミング
$ pwd   #rrbsディレクトリにいることを確認

$ trim_galore -q 20 -j 4 -rrbs rawdata/SRR1812639.fastq -o trim/
$ trim_galore -q 20 -j 4 -rrbs rawdata/SRR1812671.fastq -o trim/
$ fastqc --nogroup -o fastqc/ -t 4 trim/*.fq

$ open fastqc/SRR1812639_trimmed_fastqc.html 
$ open fastqc/SRR1812671_trimmed_fastqc.html 
  

#マッピング
$ bismark_genome_preparation --parallel 4 ~/rrbs/ref/

$ ls ref/

$ bismark ref/ trim/SRR1812639_trimmed.fq --multicore 4 -o map/
$ bismark ref/ trim/SRR1812671_trimmed.fq --multicore 4 -o map/

$ ls map/

#DNAメチル化部位の抽出

$ cd map
$ samtools sort -@ 4 SRR1812639_trimmed_bismark_bt2.bam -T tmpsam -o SRR1812639_sort.bam
$ samtools sort -@ 4 SRR1812671_trimmed_bismark_bt2.bam > SRR1812671_sort.bam

$ ls
　
$ cd ~/rrbs
$ R

#以下Rのコマンド
> library(methylKit)
> library(genomation)

> file1 <- "map/SRR1812639_sort.bam"
> file2 <- "map/SRR1812671_sort.bam"
> file.list <- list(file1,file2)

> obj <- processBismarkAln(location=file.list,sample.id=list("CLL","normal"),assembly="hg38",read.context="CpG",treatment=c(1,0))
> obj

> getMethylationStats(obj[[1]], plot=T)
> getMethylationStats(obj[[2]], plot=T)

> getCoverageStats(obj[[1]],plot=T)
> getCoverageStats(obj[[2]],plot=T)

> pdf("CLL_coverage.pdf")　   #作図デバイスを開く
> getCoverageStats(obj[[1]],plot=T)
> dev.off()                   #作図デバイスを閉じる

> meth <- unite(obj,destrand=F)
> meth 

> diff <- calculateDiffMeth(meth,mc.cores=4)
> diff

> diff25p_hyper <- getMethylDiff(diff, difference=25, qvalue=0.01, type="hyper")
> diff25p_hypo <- getMethylDiff(diff, difference=25, qvalue=0.01, type="hypo")

> gene.obj  <- readTranscriptFeatures("ref/refGene.bed")
> annotateWithGeneParts(as(diff25p_hyper,"GRanges"),gene.obj)

> annotateWithGeneParts(as(diff25p_hypo,"GRanges"),gene.obj)
　
> o1 <-order(-diff25p_hyper$meth.diff)
> o2 <-order(diff25p_hypo$meth.diff)


> a1 <- annotateWithGeneParts(as(diff25p_hyper[o1,],"GRanges"),gene.obj)
> f1 <- abs(a1@dist.to.TSS$dist.to.feature) <= 1000
> diff25p_hyper_TSS1000 <- a1@dist.to.TSS[f1,]
> head(diff25p_hyper_TSS1000) 
> write.csv(diff25p_hyper_TSS1000,"hyper.csv",quote=F)

> a2 <- annotateWithGeneParts(as(diff25p_hypo[o2,],"GRanges"),gene.obj)
> f2 <- abs(a2@dist.to.TSS$dist.to.feature) <= 1000
> diff25p_hypo_TSS1000 <- a2@dist.to.TSS[f2,]
> head(diff25p_hypo_TSS1000)

> write.csv(diff25p_hypo_TSS1000,"hypo.csv",quote=F)
> save.image("rrbs.RData")
