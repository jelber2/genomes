# Dromedary camel genome assembly

## get raw pacbio reads, (~11x haploid genome coverage)
test.sh
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR763/002/SRR7637702/SRR7637702_subreads.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR763/003/SRR7637703/SRR7637703_subreads.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR763/004/SRR7637704/SRR7637704_subreads.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR763/005/SRR7637705/SRR7637705_subreads.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR763/006/SRR7637706/SRR7637706_subreads.fastq.gz
```

```
bash test.sh > test.sh.log
```

## get error-corrected illumina reads (~65x haploid genome coverage)
test2.sh
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/003/SRR2002493/SRR2002493_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/003/SRR2002493/SRR2002493_2.fastq.gz
```
```bash
bash test2.sh
```

## error-correct pacbio reads
Snakefile-no-lighter 
```yaml
configfile: "config.yaml"
GRAPHALIGNERPATH = config["GraphAlignerPath"]
BCALMPATH = config["BcalmPath"]
BCALMCONVERTPATH = config["BcalmConvertPath"]
LIGHTERPATH = config["LighterPath"]
GENOMESIZE = config["GenomeSize"]
SHORTREADCOVERAGE = config["ShortreadCoverage"]
TMPDIR = config["TempDirectory"]
OUTDIR = config["OutputDirectory"]
SHORTREADDIR = config["ShortReadDirectory"]
SHORTREADS = config["ShortReads"]
LONGREADDIR = config["LongReadDirectory"]
LONGREADS = config["LongReads"]
SMALLK = config["SmallK"]
BIGK = config["BigK"]
ABUNDANCE = config["Abundance"]
GRAPHALIGNERPARAMS = config["GraphAlignerParams"]

SHORTREADNAMES = [n.split('.')[0] for n in SHORTREADS]
SHORTREADEXTENSIONS = ['.'.join(n.split('.')[1:]) for n in SHORTREADS]

rule all:
    input:
        OUTDIR + "corrected.fa",
        OUTDIR + "corrected_clipped.fa",
        OUTDIR + "stats.txt"

rule read_names:
    input: expand(SHORTREADDIR + "{name}.{ext}", zip, name=SHORTREADNAMES, ext=SHORTREADEXTENSIONS)
    output: temp("filenames")
    shell: "readlink -f {input} > {output}"

rule run_bcalm:
    input: 
        name = "filenames",
        files = temp(expand(TMPDIR + "{name}.cor.{ext}", zip, name=SHORTREADNAMES, ext=SHORTREADEXTENSIONS))
    output: temp("filenames.unitigs.fa")
    shadow: "full"
    log:
        stdout = TMPDIR + "bcalm_stdout.txt",
        stderr = TMPDIR + "bcalm_stderr.txt"
    threads: 40
    shell: "/usr/bin/time -v {BCALMPATH} -nb-cores {threads} -in {input.name} -kmer-size {BIGK} -abundance-min {ABUNDANCE} > {log.stdout} 2> {log.stderr}"

rule convert_bcalm:
    input: rules.run_bcalm.output
    output: TMPDIR + "graph.gfa"
    shell: "{BCALMCONVERTPATH} {input} {output} {BIGK}"

rule align_reads:
    input:
        graph = TMPDIR + "graph.gfa",
        reads = expand(LONGREADDIR + "{name}", name=LONGREADS)
    params:
        readconcat = lambda wildcards, input: ' '.join(input.reads)
    output:
        corrected = OUTDIR + "corrected.fa",
        clipped = OUTDIR + "corrected_clipped.fa"
    log:
        stdout = TMPDIR + "aligner_stdout.txt",
        stderr = TMPDIR + "aligner_stderr.txt"
    threads: 40
    shell:
        "/usr/bin/time -v {GRAPHALIGNERPATH} -g {input.graph} --corrected-out {output.corrected} --corrected-clipped-out {output.clipped} -f {params.readconcat} -t {threads} {GRAPHALIGNERPARAMS} 1> {log.stdout} 2> {log.stderr}"

rule get_stats:
    input:
        aligner_stdout = TMPDIR + "aligner_stdout.txt",
        aligner_stderr = TMPDIR + "aligner_stderr.txt",
        bcalm_stdout = TMPDIR + "bcalm_stdout.txt",
        bcalm_stderr = TMPDIR + "bcalm_stderr.txt",
        lighter_stdout = TMPDIR + "lighter_stdout.txt",
        lighter_stderr = TMPDIR + "lighter_stderr.txt"
    output:
        OUTDIR + "stats.txt"
    run:
        shell("grep 'Input reads' < {input.aligner_stdout} >> {output}")
        shell("grep 'Reads with a seed' < {input.aligner_stdout} >> {output}")
        shell("grep 'Reads with an alignment' < {input.aligner_stdout} >> {output}")
        shell("grep 'Alignments' < {input.aligner_stdout} >> {output}")
        shell("grep 'End-to-end alignments' < {input.aligner_stdout} >> {output}")
        shell("echo 'BCalm' >> {output}"),
        shell("grep 'User time' < {input.bcalm_stderr} >> {output}")
        shell("grep 'System time' < {input.bcalm_stderr} >> {output}")
        shell("grep 'Elapsed (wall clock)' < {input.bcalm_stderr} >> {output}")
        shell("grep 'Maximum resident set size' < {input.bcalm_stderr} >> {output}")
        shell("echo 'Aligner' >> {output}"),
        shell("grep 'User time' < {input.aligner_stderr} >> {output}")
        shell("grep 'System time' < {input.aligner_stderr} >> {output}")
        shell("grep 'Elapsed (wall clock)' < {input.aligner_stderr} >> {output}")
        shell("grep 'Maximum resident set size' < {input.aligner_stderr} >> {output}")
```

config.yaml
```yaml
### Change these!!
GenomeSize: 2200000000
ShortreadCoverage: 65

ShortReadDirectory: ./
# NOTE: short read endings MUST be .fq or .fa instead of .fastq or .fasta
# gzip is allowed
ShortReads:
- SRR2002493_2.fq.gz
- SRR2002493_1.fq.gz

LongReadDirectory: ./
LongReads:
- pacbio.fa

TempDirectory: tmp/
OutputDirectory: output/

# https://github.com/maickrau/GraphAligner
GraphAlignerPath: /nfs/scistore16/itgrp/jelbers/bin/GraphAligner/bin/GraphAligner
# https://github.com/GATB/bcalm
BcalmPath: /nfs/scistore16/itgrp/jelbers/bin/bcalm/build/bcalm
# https://github.com/GATB/bcalm/blob/master/scripts/convertToGFA.py
BcalmConvertPath: /nfs/scistore16/itgrp/jelbers/bin/bcalm/scripts/convertToGFA.py
# https://github.com/mourisl/Lighter
LighterPath: lighter


### Misc params. Defaults might work

# k for error correcting the reads. Try between 10-30
SmallK: 23
# k for the de Bruijn graph. Try between ~1/2 and ~2/3 of short read length
BigK: 67
# minimum k-mer abundance for the de Bruijn graph. Try between 1/100 to 2/100 of short read coverage, but not below 2.
Abundance: 3
# Parameters for GraphAligner
GraphAlignerParams: -x dbg
```

## run the correction
```bash
conda activate snakemake-7.3.6
snakemake --cores 34 all
```

## convert lowercase to uppercase bases and remove comment lines from output of GraphAligner
```bash
~/bin/seqtk/seqtk seq -UC output/corrected.fa > pacbio-GraphAligner-k67.fasta
```

## de novo assemble with flye, tried peregrine-2021, but coverage is too low
```bash
module load flye/20211026
flye --threads 96 --pacbio-corr pacbio-GraphAligner-k67.fasta --out-dir flye-pacbio-GraphAligner-k67
```

## get Hi-C reads
```bash
wget -c https://www.dropbox.com/s/k25ev41mwwxtxcs/hi-c-lib_003_R2.fq.gz
wget -c https://www.dropbox.com/s/0o6jvk66uh6oe1x/hi-c-lib_003_R1.fq.gz
wget -c https://www.dropbox.com/s/qinbgpucryu3nym/hi-c-lib_002_R2.fq.gz
wget -c https://www.dropbox.com/s/h8oa2lryrpufk88/hi-c-lib_002_R1.fq.gz
wget -c https://www.dropbox.com/s/ckrr78a4p5ywnzq/hi-c-lib_001_R2.fq.gz
wget -c https://www.dropbox.com/s/5c62ujxx5puad1q/hi-c-lib_001_R1.fq.gz
wget -c https://www.dropbox.com/s/p6eoag1k1phearn/hi-c.md5sum
md5sum --check hi-c.md5sum
```

## remove PCR duplicates
```bash
module load bcftools/1.14
module load bbtools/38.82

dedupe.sh -Xmx500g ow=t \
in1=hi-c-lib_003_R1.fq.gz usejni=t \
in2=hi-c-lib_003_R2.fq.gz \
threads=96 out=STDOUT.fa 2>dedupe.lib_003.log && \
reformat.sh int=t in=STDOUT.fa ow=t \
out1=hi-c-lib_003_R1.fa.gz \
out2=hi-c-lib_003_R2.fa.gz 2>reformat.lib_003.log &

dedupe.sh -Xmx500g ow=t \
in1=hi-c-lib_002_R1.fq.gz usejni=t \
in2=hi-c-lib_002_R2.fq.gz \
threads=96 out=STDOUT.fa 2>dedupe.lib_002.log && \
reformat.sh int=t in=STDOUT.fa ow=t \
out1=hi-c-lib_002_R1.fa.gz \
out2=hi-c-lib_002_R2.fa.gz 2>reformat.lib_002.log &

dedupe.sh -Xmx500g ow=t \
in1=hi-c-lib_001_R1.fq.gz usejni=t \
in2=hi-c-lib_001_R2.fq.gz \
threads=96 out=STDOUT.fa 2>dedupe.lib_001.log && \
reformat.sh int=t in=STDOUT.fa ow=t \
out1=hi-c-lib_001_R1.fa.gz \
out2=hi-c-lib_001_R2.fa.gz 2>reformat.lib_001.log &
```

## make symbolic links
```bash
ln -s flye-pacbio-GraphAligner-k67/assembly.fasta camel.fasta
ln -s flye-pacbio-GraphAligner-k67/assembly_graph.gfa camel.gfa
```

### error-correct camel.fasta step1
```bash
unset SLURM_EXPORT_ENV
module purge
module load bwa-mem2/2.2.1
module load samtools/1.14

bwa-mem2 mem -t 48 camel.fasta SRR2002493_1.fastq.gz SRR2002493_2.fastq.gz -R '@RG\tID:foo\tSM:bar'| \
samtools sort -@48 > SRR2002493-mapped-to-camel.fasta.bam
```

### error-correct camel.fasta step2
### call SNPs/indels with octopus (https://github.com/luntergroup/octopus)
```bash
unset SLURM_EXPORT_ENV
module purge
module load samtools/1.14
module load boost/1.78.0

samtools faidx camel.fasta
samtools index -@96 SRR2002493-mapped-to-camel.fasta.bam
/nfs/scistore16/itgrp/jelbers/bin/octopus/bin/octopus --temp temp \
-R camel.fasta -I SRR2002493-mapped-to-camel.fasta.bam \
-o SRR2002493-mapped-to-camel.fasta.bam.vcf --threads 96 --organism-ploidy 2
```

### error-correct camel.fasta step3
### use VCF file to correct errors
```bash
module load samtools/1.14
module load bcftools/1.14

bgzip -@34 SRR2002493-mapped-to-camel.fasta.bam.vcf
tabix -p vcf SRR2002493-mapped-to-camel.fasta.bam.vcf.gz

bcftools consensus --haplotype R -f camel.fasta SRR2002493-mapped-to-camel.fasta.bam.vcf.gz > camel2.fasta
```

## make sequence index for yahs (https://github.com/c-zhou/yahs) later
```bash
module load samtools/1.14
samtools faidx camel2.fasta
```

## map reads with bwa-mem2 to dustmasked genome
```bash
module load bwa-mem2/2.2.1
module load samtools/1.14
module load parallel/20220222
module load ncbi-blast/2.12.0+

function mask_data_chunk () {
  # Removes empty records and performs masking, all in pipes
  awk -v RS=">" -v FS="\n" -v ORS="" ' { if ($2) print ">"$0 } ' |\
  dustmasker -in - -outfmt fasta |\
  sed -e '/^>/!s/[a-z]/x/g'
}

export -f mask_data_chunk

cat camel2.fasta | parallel --no-notice --jobs 34 --pipe --recstart '>' \
--blocksize 100M mask_data_chunk > camel2.dustmasker.fasta

samtools faidx camel2.dustmasker.fasta

bwa-mem2 index camel2.dustmasker.fasta camel2.dustmasker.fasta > bwa-mem2-index.log 2>&1 &

bwa-mem2 mem -t 48 camel2.dustmasker.fasta <(zcat hi-c-lib_00?_R1.fa.gz) | \
perl /nfs/scistore16/itgrp/jelbers/bin/mapping_pipeline/filter_five_end.pl | \
samtools view -@48 -Sb - > hi-c1-2.bam

bwa-mem2 mem -t 48 camel2.dustmasker.fasta <(zcat hi-c-lib_00?_R2.fa.gz) | \
perl /nfs/scistore16/itgrp/jelbers/bin/mapping_pipeline/filter_five_end.pl | \
samtools view -@48 -Sb - > hi-c2-2.bam
```

## combine read 1 and read 2
```bash
# first split into 34 parts
mkdir -p partition
cd partition/
module load bbtools/38.82
partition.sh in=../hi-c1-2.bam out=hi-c1-2%.bam ways=34 2> hi-c1.partition.log &
partition.sh in=../hi-c2-2.bam out=hi-c2-2%.bam ways=34 2> hi-c2.partition.log &
```

## script for running 34 array jobs on split parts
merge.slurm

```bash
#!/bin/bash
#
#SBATCH --job-name=merge
#SBATCH -c 4
#SBATCH --mem=10g
#SBATCH --time=1:00:0
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module purge
module load samtools/1.14
module load perl/5.32.1b

i=${SLURM_ARRAY_TASK_ID}

module load samtools/1.14
perl /nfs/scistore16/itgrp/jelbers/bin/mapping_pipeline/two_read_bam_combiner.pl \
partition/hi-c1-2${i}.bam partition/hi-c2-2${i}.bam samtools 10 2>>partition/hi-c-2${i}.bam.log| \
samtools view -bS -t camel2.dustmasker.fasta.fai > partition/hi-c-2${i}.bam
```

## submit script
```bash
sbatch --array=0-33 merge.slurm
```

## script for sorting

sort.slurm

```bash
#!/bin/bash
#
#SBATCH --job-name=sort
#SBATCH -c 34
#SBATCH --mem=300g
#SBATCH --time=1:00:0
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module purge
module load samtools/1.14

samtools cat -@34 partition/hi-c-2*.bam |samtools sort -@34 -n > hi-c-2.bam
```

## submit script
```bash
sbatch sort.slurm

# if successful, then
rm -r partition/
```



## scaffold with yahs (yahs iteration 1 step)
```bash
# yahs version
cd ~/git/yahs
git show

commit e19065aa429439fc8c829687cdb5109355bca159 (HEAD -> main, tag: 1.1a.2, origin/main, origin/HEAD)

mkdir -p ~/camel/yahs
cd ~/camel/yahs/
cat ../camel2.fasta|tr '*' 'N' > camel2.fasta
module load samtools/1.14
samtools faidx camel2.fasta
~/bin/yahs/yahs -e GATC -q 10 camel2.fasta ../hi-c-2.bam > yahs.log 2>&1 &
```

## assign chromosomes
```
module load seqtk
module load ncbi-blast/2.2.31+

## first make a blastdb using BLAST 2.2.31+
makeblastdb -dbtype nucl -in yahs.out_scaffolds_final.fa

# rhmarkers/chr${i}rhmarkers.txt come from https://doi.org/10.5061/dryad.6rp36b6
mkdir rhmarkers
cd rhmarkers
wget 'https://datadryad.org/stash/downloads/file_stream/64391' -O RH-alpaca-probe-sequences.zip
unzip RH-alpaca-probe-sequences.zip
# files are chr01rhmarkers.txt,...,chr36rhmarkers.txt,chrXrhmarkers.txt (really they are FASTA files and not text files)
ls -1 chr*rhmarkers.txt|perl -pe "s/chr(\S+)rhmarkers.txt/\1/g" > samples

# change back to working directory
cd ..

# do BLAST analysis
while read i;do
## third blast the markers for each chromosome using BLAST 2.2.31+
blastn -num_threads 24 -db yahs.out_scaffolds_final.fa -query rhmarkers/chr${i}rhmarkers.txt -outfmt 6 \
-max_hsps 1 -evalue 1e-30 > rhmarkers/chr${i}rhmarkers.txt.blast
## fourth count the BLAST hits for contigs/scaffolds
echo ${i} >> contigs-to-chromosomes.txt
cut -f 2 rhmarkers/chr${i}rhmarkers.txt.blast |sort |uniq -c|sort -n |tail -n 1|perl -pe "s/( )+/\t/g" |perl -pe "s/^\t//g" |cut -f 2  >> contigs-to-chromosomes.txt
done < rhmarkers/samples &

## fifth modify numbers less than 9
perl -pi -e "s/^0//g" contigs-to-chromosomes.txt

## sixth make a copy of pilon assembly and rename it
seqtk seq -l80 yahs.out_scaffolds_final.fa > yahs.out_scaffolds_final.chromosomes.fa

##### MAKE SURE THERE ARE NO EMPTY LINES AT BEGINNING OF contigs-chromosomes.txt
## example correct format:
1
Contig0
2
Contig1

## seventh rename the contigs to chromosome names (takes about 1 hour)
cat contigs-to-chromosomes.txt | while read -r ONE;do
read -r TWO
  perl -pi -e "s/>${TWO}\n/>${ONE}\n/" yahs.out_scaffolds_final.chromosomes.fa
  echo $ONE
done

## eigth sort the chromosomes and contigs by number (ex: 1,2,3,4,X,Contig200,Contig201),
##  then output 60 bases per line, then make all bases uppercase (no soft-masking)
cat yahs.out_scaffolds_final.chromosomes.fa | seqtk seq -l0 | \
paste - - |grep -vf <(cat rhmarkers/samples|perl -pe "s/^0//g"|perl -pe "s/^/>/g") > contigs

cat yahs.out_scaffolds_final.chromosomes.fa | seqtk seq -l0 | \
paste - - |grep -f <(cat rhmarkers/samples|perl -pe "s/^0//g"|perl -pe "s/^/>/g") |fgrep -v ">X" > chromosomes

cat yahs.out_scaffolds_final.chromosomes.fa | seqtk seq -l0 | \
paste - - |fgrep ">X" > Xchromosome

cat chromosomes |sort -k 1.2 -n > tmp2 && mv tmp2 chromosomes
cat contigs |sort -k 1.8 -n > tmp2 && mv tmp2 contigs

cat chromosomes Xchromosome contigs | tr "\t" "\n" |seqtk seq -l60 -U > tmp2

mv tmp2 yahs.out_scaffolds_final.chromosomes.fa

## what percentage of the genome is not assigned to chromosomes
fgrep -v ">" yahs.out_scaffolds_final.chromosomes.fa|wc -m

2064737880

cat chromosomes Xchromosome |cut -f 2|wc -m

1753873574

cat contigs |cut -f 2|wc -m

277016954

echo "277016954/2064737880*100"|bc -l

13.41656762746077967000

# 13.4 % not assigned to chromosomes
```

## make .hic file for yahs iteration 1 step

hic.file.slurm

```bash
#!/bin/bash
#
#SBATCH --job-name=hic
#SBATCH -c 34
#SBATCH --mem=200g
#SBATCH --time=8:00:0
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module purge
module load java
module load samtools/1.14

(/nfs/scistore16/itgrp/jelbers/bin/yahs/juicer_pre yahs/yahs.out.bin yahs/yahs.out_scaffolds_final.agp yahs/camel2.fasta.fai | sort -k2,2d -k6,6d -T ./ --parallel=34 -S32G | awk 'NF' > yahs/alignments_sorted.txt.part) && (mv yahs/alignments_sorted.txt.part yahs/alignments_sorted.txt)
wget https://github.com/aidenlab/Juicebox/releases/download/v.2.13.07/juicer_tools.jar

samtools faidx yahs/yahs.out_scaffolds_final.fa
cut -f 1-2 yahs/yahs.out_scaffolds_final.fa.fai > yahs/scaffolds_final.chrom.sizes
(java -jar -Xmx200G juicer_tools.jar pre --threads 34 yahs/alignments_sorted.txt yahs/out.hic.part yahs/scaffolds_final.chrom.sizes) && (mv yahs/out.hic.part yahs/out.hic) && (rm yahs/alignments_sorted.txt)
```

## submit slurm job
```bash
sbatch hic.file.slurm
```

## generate another BAM file for yahs iteration 2
### do this in an automated way with snakePipes HiC pipeline
### https://snakepipes.readthedocs.io/en/latest/content/workflows/HiC.html#hic
### with yahs iteration 1 and hi-c reads as input
```bash
. "/nfs/scistore16/itgrp/jelbers/miniconda3/etc/profile.d/conda.sh"
conda activate snakePipes

# get FASTQ files in the working directory ~/camel/yahs
cat ../hi-c-lib_00?_R1.fq.gz > hi-c_R1.fastq.gz
cat ../hi-c-lib_00?_R2.fq.gz > hi-c_R2.fastq.gz

createIndices --local -o ./ --tools bwa \
--genomeURL /nfs/scistore16/itgrp/jelbers/camel/yahs/yahs.out_scaffolds_final.chromosomes.fa yahs.out_scaffolds_final.chromosomes > createIndices.log 2>&1

HiC \
--input-dir /nfs/scistore16/itgrp/jelbers/camel/yahs \
--output-dir /nfs/scistore16/itgrp/jelbers/camel/yahs/HiC \
--configFile /nfs/scistore16/itgrp/jelbers/miniconda3/envs/snakePipes/lib/python3.10/site-packages/snakePipes/shared/defaults.yaml \
--clusterConfigFile /nfs/scistore16/itgrp/jelbers/miniconda3/envs/snakePipes/lib/python3.10/site-packages/snakePipes/shared/cluster.yaml \
--jobs 384 \
--DAG \
--enzyme DpnII \
--noTAD \
--snakemakeOptions='--printshellcmds' \
yahs.out_scaffolds_final.chromosomes > HiC.log 2>&1
```

## get assembly quality value estimates
### use yak (https://github.com/lh3/yak)
```bash
/nfs/scistore16/itgrp/jelbers/bin/yak/yak count -b37 -t96 -o sr.yak <(zcat SRR2002493_?.fastq.gz) <(zcat SRR2002493_?.fastq.gz)
/nfs/scistore16/itgrp/jelbers/bin/yak/yak qv -t34 sr.yak camel.fasta > camel.fasta.qv.txt 2>/dev/null &

# quality values and explanation
head -n 4 camel.fasta.qv.txt && tail -n 5 camel.fasta.qv.txt 
CC	CT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count
CC	FR  fpr_lower_bound    fpr_upper_bound
CC	ER  total_input_kmers  adjusted_error_kmers
CC	CV  coverage
CT	0	0	33337680	378130.225
FR	0	0.0347
ER	2045766390	32491613.442
CV	0.989
QV	32.757	32.870

So QV is predicted between 32.757-32.870
completeness/coverage is predicted to be 98.9 %
```

## get assembly quality value estimates using the error-corrected reference
```bash
/nfs/scistore16/itgrp/jelbers/bin/yak/yak qv -t34 sr.yak camel2.fasta > camel2.fasta.qv.txt 2>/dev/null &

# quality values and explanation
head -n 4 camel2.fasta.qv.txt && tail -n 4 camel2.fasta.qv.txt 
CC	CT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count
CC	FR  fpr_lower_bound    fpr_upper_bound
CC	ER  total_input_kmers  adjusted_error_kmers
CC	CV  coverage
FR	0	0.0418
ER	2029076606	14379853.417
CV	0.993
QV	36.119	36.394

So QV is predicted between 36.119-36.394
completeness/coverage is predicted to be 99.3 %
```

## get quality values for the scaffolded assembly
```bash
/nfs/scistore16/itgrp/jelbers/bin/yak/yak qv -t34 ../sr.yak yahs.out_scaffolds_final.fa > yahs.out_scaffolds_final.qv 2> /dev/null
head -n 4 yahs.out_scaffolds_final.qv && tail -n 5 yahs.out_scaffolds_final.qv 
CC	CT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count
CC	FR  fpr_lower_bound    fpr_upper_bound
CC	ER  total_input_kmers  adjusted_error_kmers
CC	CV  coverage
CT	0	0	15314883	394767.725
FR	0	0.0418
ER	2029074776	14379824.205
CV	0.993
QV	36.119	36.394
```

## Let's view the HiC contact map made with Juicebox
the first iteration of yahs from make .hic file yahs iteration 1 step
![yahs1](https://github.com/jelber2/genomes/blob/main/pics/yahs1.svg)
