#QIIME2 workflow - determining Amplicon Sequence Variants  

**Updated 11/24/2018 - Sarah Hu**

This protocol is specific to analyzing microbial eukaryotic diversity by way of 18S rRNA gene tag sequencing. Here, we use [V4 Stoeck et al. 2010 primers](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-294X.2009.04480.x). The product is an approximately 400 bp region. [more](https://onlinelibrary.wiley.com/doi/full/10.1111/jeu.12217). Workflow below uses version 8 of QIIME 2. 

* For our QIIME1 pipeline see another [github repo](https://github.com/shu251/https://github.com/shu251/V4_tagsequencing_18Sdiversity_q1) repo.
* Molecular work related to this pipeline, for [DNA/RNA extraction protocol](https://www.protocols.io/view/rna-and-optional-dna-extraction-from-environmental-hk3b4yn)
* Also see library prep for amplifying and sequencing the V4 region of the 18S rRNA gene [here](https://www.protocols.io/view/18s-v4-tag-sequencing-pcr-amplification-and-librar-hdmb246)
* Older version of this pipeline can be found [here](https://github.com/shu251/V4_tagsequencing_18Sdiversity_qiime2)

## Requirements:
* Contents of this repo
* [QIIME2](https://docs.qiime2.org/2018.8/) version 2018.8
* R
* *See section below* Reference database to be used for downstream sequence clustering & taxonomy assignment. For 18S (microbial eukaryotic) work, I prefer [PR2](https://github.com/vaulot/pr2_database/wiki)
* To follow step by step instructions below, follow along with test files provided from here: zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1236641.svg)](https://doi.org/10.5281/zenodo.1236641)

## Activate QIIME2
First activate working environment.
```
source activate qiime2-2018.8
# Each line now starts with '(qiime2-2018.8) 
```

## Prep reference database

For several steps (including: closed OTU clusering, chimera detection, and taxonomy assignment) you will need a reference database that is imported as an artifact into qiime2.

Use below commands to set up this reference database with accompanying taxonomy information. Here, I'm using [PR2](https://github.com/vaulot/pr2_database/wiki). Download the fasta and taxonomy files suitable for mothur.

You will need to download both the fasta file and a taxonomy file. Then you'll need to import as artifacts (*create a .qza file*), [see instructions here](https://docs.qiime2.org/2018.11/tutorials/feature-classifier/)
Replace the below $PWD/pr2.fasta and pr2_tax.txt with appropriate path and reference fasta and taxonomy text file. 

You can import both the fasta DB file and taxonomy file as QIIME2 artifacts:
```
# First import the database
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $PWD/db/pr2.fasta \
  --output-path $PWD/db/pr2.qza

# Then the taxonomy file
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $PWD/db/pr2_tax.txt \
  --output-path $PWD/db/pr2_tax.qza
```
Then you have the option to select the region within the fasta db specific to the primers you used and train the classifer:
```
# Select V4 region from the PR2 database
# Use appropriate forward and reverse primers
qiime feature-classifier extract-reads \
  --i-sequences $PWD/db/pr2.qza \
  --p-f-primer CCAGCASCYGCGGTAATTCC \
  --p-r-primer ACTTTCGTTCTTGATYRA \
  --p-trunc-len 150 \
  --o-reads $PWD/db/v4_extracts.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $PWD/db/v4_extracts.qza \
  --i-reference-taxonomy $PWD/db/pr2_tax.qza \
  --o-classifier $PWD/db/pr2_classifier.qza

# tip: make sure you version the databases and taxonomy files you're using. These are often updated so you want o keep them current, but also be able to match the appropriate fasta and taxonomy file.
  
```

## Prep directories and sample IDs to import into QIIME2 environment.

This protocol assumes that you have demultiplexed paired-end data. If you have data that has not been demultiplexed, you can add a demultiplex step at the begining of the protocol, which can also be done in [QIIME2](https://docs.qiime2.org/2018.4/). QIIME2 can also handle single-end data, but the steps are not described in this protocol, as paired-end data are recommended for V4 tag sequencing analysis.

If you would like to follow along to steps below, download 'test_fastq_large.zip' files from here: zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1236641.svg)](https://doi.org/10.5281/zenodo.1236641).
These sequences are a full set of raw sequences (R1 and R2 per sample) from the SPOT station. This provides the opportunity to run the below tutorial with a full dataset.

```
unzip test_fastq_large.zip # download from zenodo
mkdir raw_seqs_dir
mv Test*.fastq* raw_seqs_dir/
```

First, create a data manifest file in a text editor, replicate the format below. Use "manifest.txt" if using the test dataset. Every line must be comma separated, with sample-id, followed by the path of fastq file, followed by either "forward" or "reverse". 

```
# Example:
sample1,$PWD/raw_seqs_dir/Test01_full_L001_R1_001.fastq.gz,forward
sample1,$PWD/raw_seqs_dir/Test01_full_L001_R2_001.fastq.gz,reverse
sample2,$PWD/raw_seqs_dir/Test02_full_L001_R1_001.fastq.gz,forward
sample2,$PWD/raw_seqs_dir/Test02_full_L001_R2_001.fastq.gz,reverse
```
* Replace $PWD with your path
* The fastq files can be gziped. 
* List all of your fastq files. 
* Save the file and name it "manifest.txt".

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33

# paired-end-demux.qza is an artifact file '.qza' that has the information stored for your raw fastq sequences

```

## Remove primers
Below script is specific for using Stoeck et al. V4 primers. Make sure to change the '--p-front-f' and '--p-front-r' input sequences if you're using other primers.  

* See [instructions](https://docs.qiime2.org/2018.8/plugins/available/cutadapt/trim-paired/)
* Use [cutadapt](https://github.com/qiime2/q2-cutadapt)
* **Citation**: Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1):pp–10, 2011. doi:10.14806/ej.17.1.200.
* Using Stoeck et al primers below. Uses IUPAC wildcards, by default cutadapt reads these
```
qiime cutadapt trim-paired \
	--i-demultiplexed-sequences paired-end-demux.qza \
	--p-cores 18 \
	--p-front-f CCAGCASCYGCGGTAATTCC \
	--p-front-r ACTTTCGTTCTTGATYRA \
	--p-error-rate 0.4 \
	--p-overlap 3 \
	--o-trimmed-sequences paired-end-demux-trimmed.qza
```

### Visualize how many reads have been removed so far:
```
# Generate the '.qzv' visual file from your '.qza' file.
## Original
qiime demux summarize --i-data paired-end-demux.qza \
	--o-visualization paired-end-demux.qzv

## After primer trimming
qiime demux summarize --i-data paired-end-demux-trimmed.qza \
	--o-visualization paired-end-demux-trimmed.qzv

# Open via local web server.
# ssh -L 8080:localhost:8080 [username]
qiime_2_ll_quick_viewer --filename paired-end-demux.qzv
qiime_2_ll_quick_viewer --filename paired-end-demux-trimmed.qzv

```

The output '.qza' (QIIME artifact) will be used directly for dada2 ASV calling (see below) or will require additional QC steps (directly below) for other OTU clustering approaches.


## Amplicon Sequence Variants (ASVs)

See [DADA2](https://docs.qiime2.org/2018.8/plugins/available/dada2/?highlight=dada2)
* Denoises PE sequences
* Dereplicates
* Filters chimera sequences out of data

```
mkdir ASVs

qiime dada2 denoise-paired \
	--i-demultiplexed-seqs paired-end-demux-trimmed.qza \
	--p-n-threads 20 \
	--p-trunc-q 2 \
	--p-trunc-len-f 200 \
	--p-trunc-len-r 200 \
	--p-max-ee 2 \
	--p-n-reads-learn 1000000 \
	--p-chimera-method pooled \
	--o-table ASVs/asv_table.qza \
	--o-representative-sequences ASVs/rep-seqs.qza \
	--o-denoising-stats ASVs/stats-dada2.qza

# Alternate options:
## Trimming, if you're having a problem with primer degradation or think cutadapt did not properly remove primers, try this instead.
	#--p-trim-left-f 20 \ 
	#--p-trim-left-r 20 \
## Chimera detection, before removing chimeras, there is an option to pool all reads ahead of detection, use the consensus where they are considered as individual sequences, or none where no chimera detection is performed.
	#--p-chimera-method
## Number of reads to learn: this is the number of reads dada2 will use to train the error model. Replace with a much smaller if you're working on testing a pipeline, this will make dada2 run much faster but it will be less accurate
```
```
#OUTPUTS:
#Saved FeatureTable[Frequency] to: ASVs/asv_table.qza
#Saved FeatureData[Sequence] to: ASVs/rep-seqs.qza
#Saved SampleData[DADA2Stats] to: ASVs/stats-dada2.qza
```

## Summarize the feature table and sequences:
```
cd ASVs

qiime feature-table summarize \
  --i-table asv_table.qza \
  --o-visualization asv_table.qzv

# In above, you can also include a metadata file, but I do not have one here ('--m-sample-metadata-file asv-metadata.tsv')
  
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

# View each of these files:
qiime_2_ll_quick_viewer --filename [xxx].qzv
```
**What to expect from visualizations:**
* *asv_table.qzv*: Give you a good breakdown of ASV output stats, i.e. how many features you have and the total frequency per sample or feature (versions of rarefaction and rank abundance plots). This interactive viewer is great in that on th second tab "Interactive Sample Detail", you can play around with sampling depth, i.e. if I assume 100,000 sequences as a base sampling depth, how many sequences do I retain? and in how many samples? The third tab, lists stats by feature, where you can see the ASV that occurred in the highest frequency, and see how many samples it appeared in. Or you can list by how many ASVs appeared in however many samples.
* *rep-seqs.qzv*: Lists the "featureID" which is a unique series of numbers and letters that is a code for what the ASV is called and then a sequence that represents that ASV, use this to download a fasta file. ALSO keep this file, you can use this to call new ASVs, but incorporate with this data. So, obviously, this would be super great to include in your supplementary data for your paper you should go write now. ALSO, you can click on each sequence and go to BLAST report page. 
* *stats-dada2.qzv*: delimited file available to download that lists each sample and how many sequences have been removed at each step, filtering, denoising, merging, and chimera-removal.

## Assign taxonomy:

Here, I am using the PR2 database - see importing and setting up reference database in the section above.
Here are two options for assigning taxonoy: to use the classifier that was trained previously or use [vsearch](https://docs.qiime2.org/2018.11/plugins/available/feature-classifier/classify-consensus-vsearch/?highlight=vsearch)

```
# Use trained classifer (See instructions above)
qiime feature-classifier classify-sklearn \
	--i-classifier /sting/tax_db/pr2_classifier.qza \
	--i-reads rep-seqs.qza \
	--o-classification asv_tax_sklearn.qza

# Or use vsearch on its own
qiime feature-classifier classify-consensus-vsearch \
	--i-query rep-seqs.qza \
	--i-reference-reads /sting/tax_db/pr2_4.10.qza \
	--i-reference-taxonomy /sting/tax_db/pr2_tax_4.10.qza \
	--p-perc-identity 0.9 \
	--p-threads 16 \
	--o-classification asv_tax_vsearch.qza
```

## Generate output tables:

```
qiime tools export \
  --input-path asv_table.qza \
  --output-path asv_table

biom convert -i feature-table.biom -o asv-table.tsv --to-tsv

qiime tools export --input-path asv_tax_vsearch.qza --output-path asv_tax_dir
# tabulated file tax_dir/taxonomy.tsv
```

Import .tsv OTU table and reformat so that column names are sample names. Compiled count information with taxonomy assignments. See script: compile_counts_tax.r


## Compile taxonomy and count tables:

Output .tsv files can be compiled in various programs. Below is an R script.

```
# Start R environment
# Import data from biom convert output
count<-read.delim("*.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

# Re-classify and re-format
colnames(count)<-as.character(unlist(count[2,]))
count1<-count[-c(1,2),]
colnames(count1)[1]<-"Feature.ID"
head(count1[1:2,]) # check dataframe structure
x<-dim(count1)[2];x # number of columns
# Convert class in count1 to numerics
count2<-count1
x<-dim(count2)[2];x # number of columns
count2[2:x] <- lapply(count2[2:x], function(x) as.numeric(as.character(x)))

# Get base stats and generate a table called "dataset_info"
seq_total<-apply(count2[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count2[2:x]>0) # OTUs per sample
OTU_single<-colSums(count2[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count2[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count2[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count2)[1] # total unique OTUs over whole dataset

# Option to remove global singletons
rowsum<-apply(count2[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1 = count2[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count2)[1] - dim(count.no1)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

# Import taxonomy file and join with count data:
taxname<-read.delim("$PWD/tax_dir/taxonomy.tsv", header=TRUE, row.names=NULL)
library(plyr) # Load plyr library into R environment
counts_wtax<-join(count.no1, taxname, by="Feature.ID", type="left", match="first")

# Write output table with OTU or ASV cluster information and taxonomy IDs:
write.table(counts_wtax, file="OutputTable_wtax.txt", quote=FALSE, sep="\t", row.names=FALSE)
```


### Get sequence run stats

In the above QC steps, whenever a .qza file was worked on, you had the option to run this:
qiime demux summarize --i-data [XXX.qza] --o-visualization XXX.qzv.  


See 'stats-dada2.qzv' to find stats of how many sequences were removed per sample:
```
qiime_2_ll_quick_viewer --filename stats-dada2.qzv
# open up in browser and see how many sequences per sample are at each step.
```


### Last updated Sarah Hu 11-24-2018

