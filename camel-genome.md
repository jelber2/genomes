# Dromedary camel genome assembly

## change directories
```bash
cd ~/camel
```

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
# but first what version of GraphAligner is being used
cd /nfs/scistore16/itgrp/jelbers/bin/GraphAligner/bin/GraphAligner
git show

commit cf7f5dbaab1e852d4a2e0c5834f56c1b49424b04 (HEAD -> master, tag: v1.0.16-osx, origin/master, origin/HEAD)
Author: Mikko Rautiainen <m_rautiainen@hotmail.com>
Date:   Fri Mar 25 15:03:28 2022 -0400

    correct mxm library in readme

cd ~/camel
conda activate snakemake-7.3.6
snakemake --snakefile Snakefile-no-lighter --cores 34 all > error-correct2.log 2>&1 &
```

## convert lowercase to uppercase bases and remove comment lines from output of GraphAligner
```bash
# https://github.com/lh3/seqtk
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

# note on octopus version
cd /nfs/scistore16/itgrp/jelbers/bin/octopus/
git show

commit 6797d8dbff78fff0900a83b2d9a3886b52602852 (HEAD -> develop, origin/develop, origin/HEAD)
Author: Daniel Cooke <daniel.cooke@invitae.com>
Date:   Fri Feb 4 12:19:06 2022 -0300

    Switch to pure Ubuntu 21.10 apt-get Docker image

```

### error-correct camel.fasta step3
### use octopus filtered VCF file to correct errors
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
### I thought dustmasking would help with repetitive regions
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

# filter_five_end.pl comes from - https://github.com/ArimaGenomics/mapping_pipeline/

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

## scaffold with yahs for iteration 2
```bash
# yahs version
cd ~/git/yahs
git show

commit e19065aa429439fc8c829687cdb5109355bca159 (HEAD -> main, tag: 1.1a.2, origin/main, origin/HEAD)

mkdir -p ~/camel/yahs2
cd ~/camel/yahs2/
ln -s ../yahs/yahs.out_scaffolds_final.chromosomes.fa yahs1.chromosomes.fasta
module load samtools/1.14
samtools faidx yahs1.chromosomes.fasta
~/bin/yahs/yahs -e GATC -q 10 yahs1.chromosomes.fasta ../yahs/HiC/HiC_matrices/hic.bam > yahs2.log 2>&1 &
```

## assign chromosomes
```
cd ~/camel/yahs2/

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

## sixth make a copy of yahs2 assembly and rename it
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
  echo $ONE >> test
done &

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

2056285562

cat chromosomes Xchromosome |cut -f 2|wc -m

1801856829

cat contigs |cut -f 2|wc -m

229043303

echo "229043303/2056285562*100"|bc -l

11.13869139737703415300

# 11.1 % not assigned to chromosomes
```

## make .hic file for yahs iteration 2 step

hic2.file.slurm

```bash
#!/bin/bash
#
#SBATCH --job-name=hic
#SBATCH -c 34
#SBATCH --mem=200g
#SBATCH --time=9:00:0
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module purge
module load java
module load samtools/1.14

(/nfs/scistore16/itgrp/jelbers/bin/yahs/juicer_pre yahs2/yahs.out.bin yahs2/yahs.out_scaffolds_final.agp yahs2/yahs1.chromosomes.fasta.fai | sort -k2,2d -k6,6d -T ./ --parallel=34 -S32G | awk 'NF' > yahs2/alignments_sorted.txt.part) && (mv yahs2/alignments_sorted.txt.part yahs2/alignments_sorted.txt)
wget -c https://github.com/aidenlab/Juicebox/releases/download/v.2.13.07/juicer_tools.jar

samtools faidx yahs2/yahs.out_scaffolds_final.fa
cut -f 1-2 yahs2/yahs.out_scaffolds_final.fa.fai > yahs2/scaffolds_final.chrom.sizes
(java -jar -Xmx200G juicer_tools.jar pre --threads 34 yahs2/alignments_sorted.txt yahs2/out3.hic.part yahs2/scaffolds_final.chrom.sizes) && (mv yahs2/out3.hic.part yahs2/out3.hic) && (rm yahs2/alignments_sorted.txt)
```

## submit slurm job
```bash
sbatch hic2.file.slurm
```

## let's trying adding the available Chicago reads for scaffolding
test3.sh

```bash
wget -c https://www.dropbox.com/s/52kq6a9od8usdbo/chicago-lib_003_R2.fq.gz
wget -c https://www.dropbox.com/s/e6lmcak6clhackc/chicago-lib_003_R1.fq.gz
wget -c https://www.dropbox.com/s/xngwermnmm4n84q/chicago-lib_002_R2.fq.gz
wget -c https://www.dropbox.com/s/1lwgr3sfbbyokcq/chicago-lib_002_R1.fq.gz
wget -c https://www.dropbox.com/s/7yjlt0lax4oq1rh/chicago-lib_001_R2.fq.gz
wget -c https://www.dropbox.com/s/vfn16yck7xohcnq/chicago-lib_001_R1.fq.gz
wget -c https://www.dropbox.com/s/hbb2wsafqyt3tlr/chicago.md5sum
md5sum --check chicago.md5sum
```

```bash
bash test3.sh > test3.sh.log 2>&1
```


## generate another BAM file for yahs iteration 3
### do this in an automated way with snakePipes HiC pipeline
### https://snakepipes.readthedocs.io/en/latest/content/workflows/HiC.html#hic
### with yahs iteration 2 and chicago and hi-c reads as input
### note that I modified the snakePipes scripts for bwa-mem2 2.2.1 instead of bwa 0.7.17
```bash
cd ~/camel/yahs2/
. "/nfs/scistore16/itgrp/jelbers/miniconda3/etc/profile.d/conda.sh"
conda activate snakePipes

# get FASTQ files in the working directory ~/camel/yahs2
cat ../chicago-lib_00?_R1.fq.gz ../hi-c-lib_00?_R1.fq.gz > chicago_hic_R1.fastq.gz
cat ../chicago-lib_00?_R2.fq.gz ../hi-c-lib_00?_R2.fq.gz > chicago_hic_R2.fastq.gz

ln -fs yahs.out_scaffolds_final.chromosomes.fa yahs2.chromosomes.fasta

createIndices --local -o ./ --tools bwa \
--genomeURL /nfs/scistore16/itgrp/jelbers/camel/yahs2/yahs2.chromosomes.fasta yahs2.chromosomes > createIndices.log 2>&1

HiC \
--input-dir /nfs/scistore16/itgrp/jelbers/camel/yahs2 \
--output-dir /nfs/scistore16/itgrp/jelbers/camel/yahs2/HiC \
--configFile /nfs/scistore16/itgrp/jelbers/miniconda3/envs/snakePipes/lib/python3.10/site-packages/snakePipes/shared/defaults.yaml \
--clusterConfigFile /nfs/scistore16/itgrp/jelbers/miniconda3/envs/snakePipes/lib/python3.10/site-packages/snakePipes/shared/cluster.yaml \
--jobs 96 \
--DAG \
--enzyme DpnII \
--noTAD \
--snakemakeOptions='--printshellcmds' \
yahs2.chromosomes > HiC.log 2>&1
```

## scaffold with yahs for iteration 3
### note: using newest yahs commit with a bug fix
```bash
# yahs version
cd ~/git/yahs
git show
commit f0803af1f347d32db79f2eb37350fc3a01506207 (HEAD -> main, origin/main, origin/HEAD)
Author: Chenxi Zhou <cz3@sanger.ac.uk>
Date:   Tue May 24 09:19:56 2022 +0100

    Bug fix

mkdir -p ~/camel/yahs3
cd ~/camel/yahs3/
ln -s ../yahs2/yahs.out_scaffolds_final.chromosomes.fa yahs2.chromosomes.fasta
module load samtools/1.14
samtools faidx yahs2.chromosomes.fasta
~/bin/yahs/yahs -e GATC -q 10 yahs2.chromosomes.fasta ../yahs2/HiC/HiC_matrices/hic.bam > yahs3.log 2>&1 &
```

## assign chromosomes
```
cd ~/camel/yahs3/

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

## sixth make a copy of yahs2 assembly and rename it
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
  echo $ONE >> test
done &

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

2064779302

cat chromosomes Xchromosome |cut -f 2|wc -m

1985599097

cat contigs |cut -f 2|wc -m

45332132

echo "45332132/2064779302*100"|bc -l

2.19549527429348475700

2.2 % not assigned to chromosomes
```

## make .hic file for yahs iteration 3 step

hic3.file.slurm

```bash
#!/bin/bash
#
#SBATCH --job-name=hic
#SBATCH -c 34
#SBATCH --mem=200g
#SBATCH --time=9:00:0
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
module purge
module load java
module load samtools/1.14

(/nfs/scistore16/itgrp/jelbers/bin/yahs/juicer_pre yahs3/yahs.out.bin yahs3/yahs.out_scaffolds_final.agp yahs3/yahs2.chromosomes.fasta.fai | sort -k2,2d -k6,6d -T ./ --parallel=34 -S32G | awk 'NF' > yahs3/alignments_sorted.txt.part) && (mv yahs3/alignments_sorted.txt.part yahs3/alignments_sorted.txt)
wget -c https://github.com/aidenlab/Juicebox/releases/download/v.2.13.07/juicer_tools.jar

samtools faidx yahs3/yahs.out_scaffolds_final.fa
cut -f 1-2 yahs3/yahs.out_scaffolds_final.fa.fai > yahs3/scaffolds_final.chrom.sizes
(java -jar -Xmx200G juicer_tools.jar pre --threads 34 yahs3/alignments_sorted.txt yahs3/out4.hic.part yahs3/scaffolds_final.chrom.sizes) && (mv yahs3/out4.hic.part yahs3/out4.hic) && (rm yahs3/alignments_sorted.txt)
```

## submit slurm job
```bash
sbatch hic3.file.slurm
```


## get assembly statistics for yahs1 and yahs2 iterations and pre-scaffolding assemblies

camel2.fasta (error-corrected pacbio reads assembled with flye --pac-corr option)

```bash
~/bin/bbmap-38.94/bbstats.sh camel2.fasta -Xmx10g
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2928	0.2076	0.2074	0.2922	0.0000	0.0000	0.0001	0.4150	0.0708

Main genome scaffold total:         	2742
Main genome contig total:           	2742
Main genome scaffold sequence total:	2030.646 MB
Main genome contig sequence total:  	2030.646 MB  	0.000% gap
Main genome scaffold N/L50:         	173/3.549 MB
Main genome contig N/L50:           	173/3.549 MB
Main genome scaffold N/L90:         	628/797.735 KB
Main genome contig N/L90:           	628/797.735 KB
Max scaffold length:                	20.748 MB
Max contig length:                  	20.748 MB
Number of scaffolds > 50 KB:        	1308
% main genome in scaffolds > 50 KB: 	99.11%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	         2,742	         2,742	 2,030,645,542	 2,030,645,542	 100.00%
    100 	         2,742	         2,742	 2,030,645,542	 2,030,645,542	 100.00%
    250 	         2,741	         2,741	 2,030,645,350	 2,030,645,350	 100.00%
    500 	         2,733	         2,733	 2,030,642,059	 2,030,642,059	 100.00%
   1 KB 	         2,611	         2,611	 2,030,555,600	 2,030,555,600	 100.00%
 2.5 KB 	         2,399	         2,399	 2,030,207,639	 2,030,207,639	 100.00%
   5 KB 	         2,197	         2,197	 2,029,463,149	 2,029,463,149	 100.00%
  10 KB 	         1,933	         1,933	 2,027,522,375	 2,027,522,375	 100.00%
  25 KB 	         1,563	         1,563	 2,021,692,341	 2,021,692,341	 100.00%
  50 KB 	         1,308	         1,308	 2,012,565,072	 2,012,565,072	 100.00%
 100 KB 	         1,118	         1,118	 1,998,995,638	 1,998,995,638	 100.00%
 250 KB 	           911	           911	 1,966,037,932	 1,966,037,932	 100.00%
 500 KB 	           761	           761	 1,912,592,259	 1,912,592,259	 100.00%
   1 MB 	           559	           559	 1,766,857,984	 1,766,857,984	 100.00%
 2.5 MB 	           266	           266	 1,289,900,333	 1,289,900,333	 100.00%
   5 MB 	            90	            90	   665,897,535	   665,897,535	 100.00%
  10 MB 	            11	            11	   149,096,641	   149,096,641	 100.00%
```

yahs1 (scaffolded camel2.fasta using yahs and Hi-C reads)

```bash
~/bin/bbmap-38.94/bbstats.sh yahs1.chromosomes.fasta -Xmx10g

A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2923	0.2075	0.2075	0.2927	0.0001	0.0000	0.0000	0.4150	0.0774

Main genome scaffold total:         	1586
Main genome contig total:           	2803
Main genome scaffold sequence total:	2030.889 MB
Main genome contig sequence total:  	2030.596 MB  	0.014% gap
Main genome scaffold N/L50:         	9/68.131 MB
Main genome contig N/L50:           	174/3.475 MB
Main genome scaffold N/L90:         	33/11.321 MB
Main genome contig N/L90:           	633/779.261 KB
Max scaffold length:                	200.47 MB
Max contig length:                  	20.748 MB
Number of scaffolds > 50 KB:        	274
% main genome in scaffolds > 50 KB: 	99.33%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	         1,586	         2,803	 2,030,888,942	 2,030,596,110	  99.99%
    100 	         1,586	         2,803	 2,030,888,942	 2,030,596,110	  99.99%
    250 	         1,585	         2,802	 2,030,888,750	 2,030,595,918	  99.99%
    500 	         1,577	         2,794	 2,030,885,459	 2,030,592,627	  99.99%
   1 KB 	         1,455	         2,672	 2,030,799,000	 2,030,506,177	  99.99%
 2.5 KB 	         1,240	         2,457	 2,030,448,039	 2,030,155,225	  99.99%
   5 KB 	         1,037	         2,254	 2,029,699,549	 2,029,406,746	  99.99%
  10 KB 	           772	         1,989	 2,027,749,775	 2,027,457,007	  99.99%
  25 KB 	           422	         1,639	 2,022,395,450	 2,022,102,754	  99.99%
  50 KB 	           274	         1,488	 2,017,184,385	 2,016,892,390	  99.99%
 100 KB 	           191	         1,387	 2,011,251,398	 2,010,963,113	  99.99%
 250 KB 	           112	         1,276	 1,998,810,889	 1,998,529,304	  99.99%
 500 KB 	            87	         1,204	 1,990,113,646	 1,989,841,721	  99.99%
   1 MB 	            69	         1,145	 1,977,919,922	 1,977,656,539	  99.99%
 2.5 MB 	            55	         1,107	 1,955,477,947	 1,955,220,025	  99.99%
   5 MB 	            44	         1,046	 1,915,227,002	 1,914,980,406	  99.99%
  10 MB 	            33	           961	 1,832,933,154	 1,832,704,075	  99.99%
  25 MB 	            21	           823	 1,616,825,718	 1,616,627,553	  99.99%
  50 MB 	            15	           734	 1,427,374,327	 1,427,197,711	  99.99%
 100 MB 	             5	           372	   731,992,790	   731,902,051	  99.99%
```

yahs2 (re-scaffolded yahs1 assembly)
```bash 
 ~/bin/bbmap-38.94/bbstats.sh yahs.out_scaffolds_final.fa -Xmx10g

A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2923	0.2075	0.2075	0.2927	0.0001	0.0000	0.0000	0.4150	0.0776

Main genome scaffold total:         	1590
Main genome contig total:           	2855
Main genome scaffold sequence total:	2030.899 MB
Main genome contig sequence total:  	2030.596 MB  	0.015% gap
Main genome scaffold N/L50:         	11/68.131 MB
Main genome contig N/L50:           	174/3.475 MB
Main genome scaffold N/L90:         	33/19.126 MB
Main genome contig N/L90:           	635/779.095 KB
Max scaffold length:                	124.144 MB
Max contig length:                  	20.748 MB
Number of scaffolds > 50 KB:        	263
% main genome in scaffolds > 50 KB: 	99.31%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	         1,590	         2,855	 2,030,898,542	 2,030,596,110	  99.99%
    100 	         1,590	         2,855	 2,030,898,542	 2,030,596,110	  99.99%
    250 	         1,589	         2,854	 2,030,898,350	 2,030,595,918	  99.99%
    500 	         1,581	         2,846	 2,030,895,059	 2,030,592,627	  99.99%
   1 KB 	         1,459	         2,724	 2,030,808,600	 2,030,506,177	  99.99%
 2.5 KB 	         1,243	         2,508	 2,030,456,639	 2,030,154,225	  99.99%
   5 KB 	         1,040	         2,305	 2,029,708,149	 2,029,405,746	  99.99%
  10 KB 	           774	         2,039	 2,027,750,375	 2,027,448,007	  99.99%
  25 KB 	           417	         1,681	 2,022,255,152	 2,021,953,068	  99.99%
  50 KB 	           263	         1,524	 2,016,878,244	 2,016,576,868	  99.99%
 100 KB 	           163	         1,404	 2,009,732,182	 2,009,434,931	  99.99%
 250 KB 	            98	         1,315	 1,999,689,854	 1,999,397,644	  99.99%
 500 KB 	            77	         1,255	 1,992,140,197	 1,991,856,022	  99.99%
   1 MB 	            64	         1,215	 1,982,627,755	 1,982,349,228	  99.99%
 2.5 MB 	            52	         1,181	 1,962,849,063	 1,962,575,523	  99.99%
   5 MB 	            44	         1,139	 1,932,638,975	 1,932,373,397	  99.99%
  10 MB 	            37	         1,085	 1,879,539,021	 1,879,284,468	  99.99%
  25 MB 	            25	           911	 1,652,712,350	 1,652,496,667	  99.99%
  50 MB 	            17	           757	 1,385,733,105	 1,385,553,385	  99.99%
 100 MB 	             4	           219	   463,877,766	   463,824,343	  99.99%
```

yahs3 (re-scaffolded yahs2 assembly)

```bash
~/bin/bbmap-38.94/bbstats.sh yahs.out_scaffolds_final.fa
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2924	0.2074	0.2076	0.2926	0.0002	0.0000	0.0000	0.4150	0.0792

Main genome scaffold total:         	1487
Main genome contig total:           	2908
Main genome scaffold sequence total:	2030.930 MB
Main genome contig sequence total:  	2030.596 MB  	0.016% gap
Main genome scaffold N/L50:         	11/74.119 MB
Main genome contig N/L50:           	176/3.428 MB
Main genome scaffold N/L90:         	29/24.248 MB
Main genome contig N/L90:           	654/740.176 KB
Max scaffold length:                	124.144 MB
Max contig length:                  	20.748 MB
Number of scaffolds > 50 KB:        	183
% main genome in scaffolds > 50 KB: 	99.35%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	         1,487	         2,908	 2,030,929,742	 2,030,596,110	  99.98%
    100 	         1,487	         2,908	 2,030,929,742	 2,030,596,110	  99.98%
    250 	         1,486	         2,907	 2,030,929,550	 2,030,595,918	  99.98%
    500 	         1,478	         2,899	 2,030,926,259	 2,030,592,627	  99.98%
   1 KB 	         1,356	         2,777	 2,030,839,800	 2,030,506,177	  99.98%
 2.5 KB 	         1,140	         2,561	 2,030,487,839	 2,030,154,225	  99.98%
   5 KB 	           936	         2,357	 2,029,735,349	 2,029,401,746	  99.98%
  10 KB 	           670	         2,091	 2,027,777,575	 2,027,444,007	  99.98%
  25 KB 	           316	         1,736	 2,022,351,459	 2,022,018,178	  99.98%
  50 KB 	           183	         1,598	 2,017,761,322	 2,017,429,098	  99.98%
 100 KB 	           110	         1,507	 2,012,497,537	 2,012,168,997	  99.98%
 250 KB 	            64	         1,425	 2,005,173,176	 2,004,852,025	  99.98%
 500 KB 	            52	         1,385	 2,001,004,807	 2,000,689,367	  99.98%
   1 MB 	            44	         1,345	 1,995,643,938	 1,995,335,043	  99.98%
 2.5 MB 	            37	         1,280	 1,985,599,060	 1,985,302,033	  99.99%
   5 MB 	            36	         1,269	 1,981,241,365	 1,980,946,622	  99.99%
  10 MB 	            36	         1,269	 1,981,241,365	 1,980,946,622	  99.99%
  25 MB 	            28	         1,141	 1,810,708,287	 1,810,442,780	  99.99%
  50 MB 	            17	           901	 1,435,381,973	 1,435,172,395	  99.99%
 100 MB 	             4	           308	   468,884,874	   468,814,310	  99.98%
```

## get assembly quality value estimates
### use yak (https://github.com/lh3/yak)
first assembly without error-correction by correcting SNPs/indels called with octopus
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
second assembly with SNPs/indels corrected by bcftools consensus based on octopus called variants
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
#### results are pretty much identical between yahs1, yahs2, and yahs3 iterations
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
the first iteration of yahs from [make .hic file yahs iteration 1 step](https://github.com/jelber2/genomes/blob/main/camel-genome.md#scaffold-with-yahs-yahs-iteration-1-step)
[yahs1](https://github.com/jelber2/genomes/blob/main/pics/yahs1.svg)

the second iteration of yahs from [make .hic file yahs iteration 2 step](https://github.com/jelber2/genomes/blob/main/camel-genome.md#make-hic-file-for-yahs-iteration-2-step)
[yahs2](https://github.com/jelber2/genomes/blob/main/pics/yahs2.svg)

the third iteration of yahs from [make .hic file yahs iteration 3 step](https://github.com/jelber2/genomes/blob/main/camel-genome.md#make-hic-file-for-yahs-iteration-3-step)
[yahs3](https://github.com/jelber2/genomes/blob/main/pics/yahs3.svg)


## Let's view a dotplot of a whole genome alignment between CamDro3(https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000803125.2/) and yahs3.chromosomes.fasta

[dotplot](https://github.com/jelber2/genomes/blob/main/pics/map_yahs3_to_CamDro3.svg)
