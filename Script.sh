

#Set working environment 

#Install bioconda 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  ## follow the prompt. Keep pressing ENTER or responding by yes when needed :)
## restart the terminal
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda

#Create conda env
conda create -y --name ngs1 python=3.6

################################################### Install the software #####################################################################
conda activate ngs1
#For Quality Control
conda install -c bioconda fastqc 
conda install -c bioconda multiqc 
#For Alignmnent 
conda install -c bioconda -y bwa
# install Samtools
conda install -y samtools
# Install Picard tools & make sure you have the right Java version
conda install -c bioconda picard 
picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0
conda install -c bioconda java-jdk=8.0.112
java -version
#Install GATK
conda install -c bioconda gatk4 
## install ggplot2 required for the AnalyzeCovariates tool to plot the QC plots 
conda install -c r r-ggplot2
conda install -c conda-forge r-gsalib
#VCFtool
 conda install -c bioconda vcftools 





# Ceate work directory 

mkdir ~/Breast_cancer_samples && cd ~/Breast_cancer_samples
############################################################# DATA ###################################################################################
#Download Samples (each sample is 2 reads ) 
 curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR730/002/SRR7309332/SRR7309332_1.fastq.gz -o SRR7309332_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_1.fastq.gz
 curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR730/002/SRR7309332/SRR7309332_2.fastq.gz -o SRR7309332_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_2.fastq.gz
 curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR730/008/SRR7309338/SRR7309338_1.fastq.gz -o SRR7309338_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_1.fastq.gz
 curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR730/008/SRR7309338/SRR7309338_2.fastq.gz -o SRR7309338_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_2.fastq.gz
 curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR730/005/SRR7309325/SRR7309325_1.fastq.gz -o SRR7309325_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_1.fastq.gz
 curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR730/005/SRR7309325/SRR7309325_2.fastq.gz -o SRR7309325_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_2.fastq.gz

#Download reference:
 
''' download chr2, 7, 12, 17 after that concatenate them in one fasta file (https://www.biostars.org/p/256796/) (https://unix.stackexchange.com/questions/158941/how-to-combine-gunzipped-fastq-files) '''


wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.12.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.chromosome.12.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa.gz
cat Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa Homo_sapiens.GRCh38.dna_sm.chromosome.7.fa Homo_sapiens.GRCh38.dna_sm.chromosome.12.fa Homo_sapiens.GRCh38.dna_sm.chromosome.17.fa > Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa 

########################################################## Methods ################################################################# 
#############################
#check the quality of data
#Run the FASTQC for each read end
mkdir ~/Breast_cancer_samples/FASTQC_sample1_2 && cd ~/Breast_cancer_samples/FASTQC_sample1_2
for f in ~/Breast_cancer_samples/*.fastq.gz;do fastqc -t 1 -f fastq -noextract $f;done
#Merge the output reports into one super report
mv ../fqData/*html ./
mv ../fqData/*zip ./
multiqc -z -o . .
############################ 
#alligment 

#index your genome
mkdir -p  workdir/bwa_align/bwaIndex && cd workdir/bwa_align/bwaIndex

ln -s ~/Desktop/Breast_cancer_samples/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa .
bwa index -a bwtsw Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa

mkdir -p ~/Desktop/Breast_cancer_samples/GATK_varient_calling && cd ~/Desktop/Breast_cancer_samples/GATK_varient_calling

# renam the sample name to be seutable to gtak 

mv SRR7309325_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_1.fastq.gz SRR7309325_L001_R1.fastq.gz
mv SRR7309325_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_2.fastq.gz SRR7309325_L001_R2.fastq.gz
mv SRR7309332_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_1.fastq.gz SRR7309332_L002_R1.fastq.gz  
mv SRR7309332_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_2.fastq.gz SRR7309332_L002_R2.fastq.gz    
mv SRR7309338_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_1.fastq.gz SRR7309338_L002_R1.fastq.gz
mv SRR7309338_DNA-seq_of_Adult_Female_HER2_Breast_Cancer_2.fastq.gz SRR7309338_L002_R2.fastq.gz

for R1 in ~/Desktop/Breast_cancer_samples/*_R1.fastq.gz;do
    SM=$(basename $R1 | cut -d"_" -f1)                                          ##sample ID
    LB=$(basename $R1 | cut -d"_" -f1,2)                                        ##library ID
    PL="Illumina"                                                               ##platform (e.g. illumina, solid)
    RGID=$(zcat $R1 | head -n1 | sed 's/ /_/g' |cut -d "_" -f1)                 ##read group identifier 
    PU=$RGID.$LB                                                                ##Platform Unit
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"
     R2=$(echo $R1 | sed 's/_R1./_R2./')
    echo $R1 $R2
    index="/home/admin1/Desktop/Breast_cancer_samples/workdir/bwa_align/bwaIndex/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa"
    bwa mem -t 4 -M -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" $index $R1 $R2 > $(basename $R1 _R1_001.pe.fq.gz).sam
done

#see statistics of the alligment 

samtools flagstat ~/Desktop/Breast_cancer_samples/GATK_varient_calling/SRR7309325_L001_R1.fastq.gz.sam > SRR7309325_L001_R1.fastq.gz.sam_stats.out
samtools flagstat ~/Desktop/Breast_cancer_samples/GATK_varient_calling/SRR7309332_L002_R1.fastq.gz.sam > SRR7309332_L002_R1.fastq.gz.sam_stats.out
samtools flagstat ~/Desktop/Breast_cancer_samples/GATK_varient_calling/SRR7309338_L002_R1.fastq.gz.sam > SRR7309338_L002_R1.fastq.gz.sam_stats.out


#generate & sort BAM file

samtools view -hbo SRR7309325_L001_R1.fastq.gz.bam SRR7309325_L001_R1.fastq.gz.sam
samtools sort SRR7309325_L001_R1.fastq.gz.bam -o SRR7309325_L001_R1.fastq.gz.sorted.bam

samtools view -hbo SRR7309332_L002_R1.fastq.gz.bam SRR7309332_L002_R1.fastq.gz.sam
samtools sort SRR7309332_L002_R1.fastq.gz.bam -o SRR7309332_L002_R1.fastq.gz.sorted.bam


samtools view -hbo SRR7309338_L002_R1.fastq.gz.bam SRR7309338_L002_R1.fastq.gz.sam
samtools sort SRR7309338_L002_R1.fastq.gz.bam -o SRR7309338_L002_R1.fastq.gz.sorted.bam



#############################
#mark deduplicates

picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0
 
for sample in *.sorted.bam;do
  name=${sample%.sorted.bam}
  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done

samtools flagstat SRR7309325_L001_R1.fastq.gz.dedup.bam > SRR7309325_L001_R1.fastq.gz.dedup.stat
samtools flagstat SRR7309332_L002_R1.fastq.gz.dedup.bam > SRR7309332_L002_R1.fastq.gz.dedup.stat
samtools flagstat SRR7309338_L002_R1.fastq.gz.dedup.bam > SRR7309338_L002_R1.fastq.gz.dedup.stat

##############################
#indexing to gtak variant calling 

#a) samples

java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR7309325_L001_R1.fastq.gz.dedup.bam
java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR7309332_L002_R1.fastq.gz.dedup.bam
java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=SRR7309338_L002_R1.fastq.gz.dedup.bam

#b) reference 

ln -s /home/admin1/Desktop/Breast_cancer_samples/Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa O=Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.dict
samtools faidx Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa

#Download known varinats & concatenat them in  a one file and index this file: 

wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr2.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr7.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr12.vcf.gz
wget -c ftp://ftp.ensembl.org/pub/release-99/variation/vcf/homo_sapiens/homo_sapiens-chr17.vcf.gz
gunzip  homo_sapiens-chr2.vcf.gz
gunzip  homo_sapiens-chr7.vcf.gz
gunzip  homo_sapiens-chr12.vcf.gz
gunzip  homo_sapiens-chr17.vcf.gz
#concatenate all on one file  
grep "^#" homo_sapiens-chr2.vcf > homo_sapiens-chr2_7_12_17.vcf
grep "^2" homo_sapiens-chr2.vcf >> homo_sapiens-chr2_7_12_17.vcf
grep "^7" homo_sapiens-chr7.vcf >> homo_sapiens-chr2_7_12_17.vcf
grep "^12" homo_sapiens-chr12.vcf >> homo_sapiens-chr2_7_12_17.vcf
grep "^17" homo_sapiens-chr17.vcf >> homo_sapiens-chr2_7_12_17.vcf
gatk IndexFeatureFile -I homo_sapiens-chr2_7_12_17.vcf

##############################
#recalibration 
#Recalibrate Bases BQSR

 
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa  -I $sample --known-sites homo_sapiens-chr2_7_12_17.vcf \
-O $name.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I $name.bqsr.bam --known-sites homo_sapiens-chr2_7_12_17.vcf \
-O $name.report2

  gatk AnalyzeCovariates \
-before $name.report \
-after $name.report2 \
-plots $name.pdf
done

############################## ?????

#Joint variant calling using HaplotypeCaller

#Call germline SNPs and indels via local re-assembly of haplotypes

# assess genotype likelihood per-sample
for sample in *.bqsr.bam;do
  name=${sample%.bqsr.bam}

  gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -I $sample \
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -O $name.gvcf
done

## combine samples
gatk --java-options "-Xmx2G" CombineGVCFs -R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa -V SRR7309325_L001_R1.fastq.gz.gvcf -V SRR7309332_L002_R1.fastq.gz.gvcf  \
-V SRR7309338_L002_R1.fastq.gz.gvcf  -O raw_variants.gvcf

## Joint Genotyping
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 6 -O raw_variants.vcf

## annotated output
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
--dbsnp homo_sapiens-chr2_7_12_17.vcf \
-O raw_variants_ann.vcf

## check how many variant got annotated
grep -v "^#" raw_variants_ann.vcf | awk '{print $3}' | grep "^rs" | wc -l

##############################
#VCF statitics
#First letus index the VCF file
conda install -c bioconda tabix
bgzip -c raw_variants_ann.vcf > raw_variants_ann.vcf.gz
tabix -p vcf raw_variants_ann.vcf.gz
#Calc some stats about your vcf
conda install -c bioconda rtg-tools
rtg vcfstats raw_variants_ann.vcf.gz > stats.txt

##############################
#Split SNPs and indels

gatk --java-options "-Xmx2G" SelectVariants \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants_ann.vcf \
--select-type-to-include SNP \
-O raw_variants_ann_SNP.vcf

gatk --java-options "-Xmx2G" SelectVariants \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants_ann.vcf \
--select-type-to-include INDEL \
-O raw_variants_ann_INDEL.vcf

##############################

#SNP Variant filteration


gatk --java-options "-Xmx2G" VariantFiltration \
-R Homo_sapiens.GRCh38.dna_sm.chromosome2_7_12_17.fa \
-V raw_variants_ann_SNP.vcf \
--filter-name "snpQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "snpMQ" \
--filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filter-name "snpMQRankSum" \
--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filter-name "snpFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
--filter-name "snpSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filter-name "snpReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filter-name "snpDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_SNP_clean.vcf




# Assess the different filters in both known and novel


for var in "SNP" "INDEL";do
 input="raw_variants_ann_"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done;

mkdir filters && cd filters
mv ../{*.SNP.*,SNP.*,*.INDEL.*,INDEL.*} .

wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done
done 


# extract all passed record in a new vcf
vcftools --vcf raw_variants_ann_SNP_clean.vcf --remove-filtered-all --recode --recode-INFO-all --stdout > raw_variants_ann_SNP_clean_all_passed.vcf

#VCF statitics to the passed file

bgzip -c raw_variants_ann_SNP_clean_all_passed.vcf > raw_variants_ann_SNP_clean_all_passed.vcf.gz
tabix -p vcf raw_variants_ann_SNP_clean_all_passed.vcf.gz

rtg vcfstats raw_variants_ann_SNP_clean_all_passed.vcf.gz > raw_variants_ann_SNP_clean_all_passed.stats.txt

#extract known snps for further analysis in a text file:

grep -v "^#" raw_variants_ann_SNP_clean_all_passed.vcf | awk '{print $3}' | grep "^rs" > raw_variants_ann_SNP_clean_all_passed.Id.txt
