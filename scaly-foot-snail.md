# Analysis of Chrysomallon squamiferum

## Main directory

```sh
cd ~/test/
```

## Get original Illumina WGS reads
```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/007/SRR8599727/SRR8599727_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/007/SRR8599727/SRR8599727_2.fastq.gz
```

## Get Hi-C reads
```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/009/SRR8599719/SRR8599719_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR859/009/SRR8599719/SRR8599719_2.fastq.gz
```

## Get rebase-called Nanopore reads
```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/091/SRR12763791/SRR12763791_1.fastq.gz
```

## Decontaminate Illumina reads

decontaminate-snakefile

```sh
#
shell.executable("/bin/bash")

# set number of threads here
THREADS=48
IDS, = glob_wildcards("{id}_2.fastq.gz")

# list of rules which are not deployed to slurm
localrules: all, taxids, getLinks, download

rule all:
       input: expand("{id}_decon_1.fastq.gz", id = IDS)

# Use sendsketch.sh 38.82 to assess potential contamination
# from BBMap/BBTools https://sourceforge.net/projects/bbmap/ .
#
# Note that an alternative to the `sendsketch.sh`
# steps is to download
# RefSeq bacteria, RefSeq viruses, RefSeq Human
# for example then run `seal.sh`
# to remove those sequences from the HiFi reads
rule sendsketch:
        input: 
            read1=expand("{id}_1.fastq.gz", id = IDS),
            read2=expand("{id}_2.fastq.gz", id = IDS)
        output: 
            sketch="sketch.sketch",
            log="sendsketch.log.txt"
        shell:
            """
            module purge
            module load compression-tools/20220329
            module load bbtools/38.82
            reformat.sh in1={input.read1} in2={input.read2} out=STDOUT.fa | \
            sendsketch.sh in=STDIN.fa threads={THREADS} \
            outsketch={output.sketch} \
            out={output.log} records=1000
            """

# get genome taxids with genome size < 35 Mbp that appear as contaminants
rule taxids:
        input: "sendsketch.log.txt"
        output: "genomes/genomes.to.get"
        shell:
            """
            tail -n +4 {input} | \
            cut -f 9,10,11,12| \
            perl -pe 's/K(\s+)/000$1/g' | \
            perl -pe 's/M(\s+)/000000$1/g'| \
            awk '$2 < 35000000' OFS='\t' |cut -f 1 | \
            sort -u |grep -P '^\d' > {output}
            """

# use NCBI Datasets to retrieve genomes from taxids
# get the links in a zip file for each genome
rule getLinks:
        input: "genomes/genomes.to.get"
        output: "genomes/to.download.txt2"
        params: "genomes"
        shell:
            """
            # get dehydrated files
            wget -c https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
            chmod u+x datasets
            while read i
            do
              ./datasets download genome taxon $i --exclude-gff3 --reference \
              --exclude-protein --exclude-rna \
              --dehydrated --filename {params}/$i.zip || :
            done < {input}

            # get links
            for i in `ls {params}/*.zip`;do
              unzip -p $i *fetch.txt | \
              grep 'GCF' |head -n 1|grep -v 'json' >> {params}/to.download.txt || :
            done

            # get links to genome
            cut -f 1 {params}/to.download.txt > {output}
            """

#### actually download the FASTA files into the file genomes.fasta.gz
rule download:
        input: "genomes/to.download.txt2"
        output: "genomes/genomes.fasta.gz"
        params: "genomes"
        shell:
            """
            while read i
            do
              wget $i -O - >> {output} || :
            done < {input}
            """

### Use seal.sh 38.82
#from BBMap/BBTools https://sourceforge.net/projects/bbmap/ .
#Use seal.sh with minkmerfraction=0.5 (50% 31-mers from reads must match 
#references to be contamination using libraries of "contamination"
#based on sendsketch.sh). seal.sh seems to have better accuracy than kraken2

rule seal:
        input:
            read1="{id}_1.fastq.gz",
            read2="{id}_2.fastq.gz",
            ref="genomes/genomes.fasta.gz"
        output:
            read1="{id}_decon_1.fastq.gz",
            read2="{id}_decon_2.fastq.gz"
        shell:
            """
            module purge
            module load compression-tools/20220329
            module load bbtools/38.82
            seal.sh threads={THREADS} k=31 ow=t \
            ref=genomes/genomes.fasta.gz \
            minkmerfraction=0.5 \
            in1={input.read1} in2={input.read2} \
            outu1={output.read1} outu2={output.read2} \
            out=/dev/null
            """
```

perform decontamination

```sh
snakemake -j 10 --snakefile decontaminate-snakefile --printshellcmds --latency-wait 60 --local-cores 4 \
--cores all --cluster "sbatch --export=NONE --no-requeue --job-name {rule} --mem=100g \
--time=4:00:00 --cpus-per-task=48 " all > decontaminate-snakefile.log 2>&1 &
```


### Results of decontamination

```sh
Input:                          356153620 reads    53423043000 bases.
Matched reads:                  16386480  reads    2457972000 bases (4.60%)
Unmatched reads:                339767140 reads   50965071000 bases (95.40%)
```


examples of microrganism contamination (not just bacteria)

```sh
# Genus
Achromobacter
Anaerostipes
Ascoidea
Aspergillus
Bacillus
Bacteroides
Blautia
Cavenderia
Chlamydia
Clostridium
Dictyostelium
Entamoeba
Enterobacter
Escherichia
Eubacterium
Holdemanella
Mycobacterium
Mycobacteroides
Mycolicibacterium
Naegleria
Plasmodium
Ruminococcus
Shigella
```


## Subsample Nanopore reads to 50x assuming 400M genome

use https://github.com/mbhall88/rasusa version 0.6.1

```sh
~/bin/rasusa-0.6.1-x86_64-unknown-linux-musl/rasusa --seed 1 --coverage 50 \
--genome-size 400M --input SRR12763791_1.fastq.gz \
2> rasusa.log |seqtk seq -A > SRR12763791.50x.fasta
```

## Quality and adapter trimming
```sh
module load bbtools/38.82
bbduk.sh threads=24 in1=SRR8599727_decon_1.fastq.gz in2=SRR8599727_decon_2.fastq.gz \
out1=SRR8599727_decon_trim_1.fastq.gz out2=SRR8599727_decon_trim_2.fastq.gz \
ref=~/bin/bbmap-38.94/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15
```

## Error-correct decontaminated, QC'd Illumina reads
```sh
module load bcftools/1.14
module load bbtools/38.82

error-correct using k=31 and tadpole.sh 

tadpole.sh threads=34 in1=SRR8599727_decon_trim_1.fastq.gz in2=SRR8599727_decon_trim_2.fastq.gz \
out1=SRR8599727_decon_trim_corr_1.fastq.gz out2=SRR8599727_decon_trim_corr_2.fastq.gz mode=correct
```

## Interleave for KmerGenie
http://kmergenie.bx.psu.edu/

```sh
module load bbtools/38.82
reformat.sh in=SRR8599727_decon_trim_corr_1.fastq.gz in2=SRR8599727_decon_trim_corr_2.fastq.gz \
out=SRR8599727_decon_trim_corr_interleaved.fasta

/nfs/scistore16/itgrp/jelbers/bin/kmergenie-1.7051/kmergenie \
SRR8599727_decon_trim_corr_interleaved.fasta -t 48 -k 121 -l 61
```

'best' K=111
![KmerGenie](https://github.com/jelber2/genomes/blob/main/pics/KmerGenie.png)


## Use GraphAligner commit # cf7f5db and decontaminated, QC'd, and error-corrected Illumina WGS reads to correct 50x Nanopore

first make unitigs with BCALM 2, version v2.2.3, git commit b8cde9c (https://github.com/GATB/bcalm)

```sh
ls -1 SRR8599727_decon_trim_corr_1.fastq.gz SRR8599727_decon_trim_corr_2.fastq.gz > filenames
/nfs/scistore16/itgrp/jelbers/bin/bcalm/build/bcalm -nb-cores 48 -in filenames -kmer-size 111 -abundance-min 3
```

then convert unitigs.fa to untigs.gfa

```sh
/nfs/scistore16/itgrp/jelbers/bin/bcalm/scripts/convertToGFA.py filenames.unitigs.fa filenames.unitigs.gfa 111
```

now error-correct the long reads

```sh
/nfs/scistore16/itgrp/jelbers/bin/GraphAligner/bin/GraphAligner -g filenames.unitigs.gfa \
--corrected-out corrected.50x.fa -f SRR12763791.50x.fasta -t 96 -x dbg
```
  
make lower case bases upper case and put on a single line

```sh
# seqtk 1.3-r117-dirty
seqtk seq -Ul0 corrected.50x.fa > SRR12763791.50x.GraphAligner2.fasta
```

## Assemble with flye (2.9-b1778)

```sh
/nfs/scistore16/itgrp/jelbers/git/Flye/bin/flye --threads 48 \
--nano-corr SRR12763791.50x.GraphAligner2.fasta --out-dir flye-SRR12763791.50x.GraphAligner
```

## Purge dups commit # 44fcd38

```sh
cd flye-SRR12763791.50x.GraphAligner

export "PATH=/nfs/scistore16/itgrp/jelbers/bin/purge_dups/src:$PATH"
split_fa p_ctgs.fasta > split.fa

# mm2-fast commit # 830e8c7 from https://github.com/bwa-mem2/mm2-fast
/nfs/scistore16/itgrp/jelbers/bin/mm2-fast/minimap2 -t 24 -x map-ont --secondary=no \
--max-chain-skip=1000000 assembly.fasta ../SRR12763791.50x.GraphAligner2.fasta \
2>/dev/null | pigz -p 24 -c - > test.paf.gz

/nfs/scistore16/itgrp/jelbers/bin/mm2-fast/minimap2 --max-chain-skip=1000000 \
-t 24 -xasm5 -DP split.fa split.fa 2>/dev/null| pigz -p 24 -c - > split.fa.paf.gz

pbcstat test.paf.gz -O ./
calcuts PB.stat > PB.cutoffs
purge_dups -2 -T PB.cutoffs -c PB.base.cov split.fa.paf.gz > PB.bed
get_seqs -p assembly PB.bed assembly.fasta
```

## Hi-C

quality and adapter trimming of Hi-C reads

```sh
module load bbtools/38.82
bbduk.sh threads=24 in1=SRR8599719_1.fastq.gz in2=SRR8599719_2.fastq.gz \
out1=SRR8599719_trim_1.fastq.gz out2=SRR8599719_trim_2.fastq.gz \
ref=~/bin/bbmap-38.94/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 \
tpe tbo qtrim=rl trimq=15 > bbduk.log 2>&1
```

remove PCR duplicates

```sh
module load bbtools/38.82
dedupe.sh usejni=t in=SRR8599719_trim_1.fastq.gz in2=SRR8599719_trim_2.fastq.gz \
out=SRR8599719_trim_dedup_interleaved.fastq.gz threads=48

# interleaved to pairs
module load bbtools/38.82
reformat.sh int=t in=SRR8599719_trim_dedup_interleaved.fastq.gz out1=SRR8599719_trim_dedup_1.fastq.gz \
out2=SRR8599719_trim_dedup_2.fastq.gz
```

map reads with BWA-MEM2, filter and combine with https://github.com/ArimaGenomics/mapping_pipeline/

```sh
module load bwa-mem2/2.2.1
module load samtools/1.14
module load perl/5.32.1b
cd ~/test
cp flye-SRR12763791.50x.GraphAligner/assembly.purged.fa assembly.purged.fasta
bwa-mem2 index -p assembly.purged.fasta assembly.purged.fasta 2>/dev/null &

module purge
module load bwa-mem2/2.2.1
module load samtools/1.14
module load perl/5.32.1b
bwa-mem2 mem -t 96 assembly.purged.fasta SRR8599719_trim_dedup_1.fastq.gz | \
perl /nfs/scistore16/itgrp/jelbers/bin/mapping_pipeline/filter_five_end.pl | \
samtools view -Sb > R1.bam

module purge
module load bwa-mem2/2.2.1
module load samtools/1.14
module load perl/5.32.1b
bwa-mem2 mem -t 96 assembly.purged.fasta SRR8599719_trim_dedup_2.fastq.gz | \
perl /nfs/scistore16/itgrp/jelbers/bin/mapping_pipeline/filter_five_end.pl | \
samtools view -Sb > R2.bam
```

keep only primary alignments

```sh
module load samtools/1.14
samtools view -h -@12 -F 2304 R1.bam |samtools view -@8 -Sb > R2b.bam  &
samtools view -h -@12 -F 2304 R2.bam |samtools view -@8 -Sb > R2b.bam  &
```

split into 64 pieces so we can do parallel processing to combine R1 and R2

```sh
cd ~/test
mkdir -p partition
cd partition
partition.sh in=../R1b.bam out=hi-c1-2%.bam ways=64 2> hi-c1.partition.log &
partition.sh in=../R2b.bam out=hi-c2-2%.bam ways=64 2> hi-c2.partition.log &
```

### script for parallel processing

merge.slurm

```bash
#!/bin/bash
#
#SBATCH --job-name=merge
#SBATCH -c 4
#SBATCH --mem=12g
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
partition/hi-c1-2${i}.bam partition/hi-c2-2${i}.bam samtools 15 2>>partition/hi-c-2${i}.bam.log| \
samtools sort > partition/hi-c-2${i}.bam
```

submit array job

```sh
sbatch --array=0-63 merge.slurm
```

sort by read name

```sh
module load samtools/1.14
samtools cat -@34 partition/hi-c-2*.bam |samtools sort -@34 -n > hi-c-2.bam
```

## Scaffold

https://github.com/c-zhou/yahs
YaHS commit# f0803af

```sh
module load samtools/1.14

cd ~/test
mkdir -p assembly.purged
cd assembly.purged
cp ../assembly.purged.fasta assembly.purged.fasta
samtools faidx assembly.purged.fasta
~/bin/yahs/yahs -e AAGCTT assembly.purged.fasta hi-c-2.bam > yahs.log 2>&1 &
```

make Hi-C contact map

```sh
cd ~/test/assembly.purged

(/nfs/scistore16/itgrp/jelbers/bin/yahs/juicer_pre yahs.out.bin yahs.out_scaffolds_final.agp assembly.purged.fasta.fai | sort -k2,2d -k6,6d -T ./ --parallel=24 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)
wget -c https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar

samtools faidx yahs.out_scaffolds_final.fa
cut -f 1-2 yahs.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes
(java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre --threads 24 alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) && (mv out4.hic.part out4.hic)
```

Dot plot

```sh
$HOME/bin/mm2-fast/minimap2 --max-chain-skip=1000000 -cx asm5 -t24 --cs ~/test/Dryad_upload/Csq_v2.0.fasta.gz yahs.out_scaffolds_final.chromosomes.fa 2>/dev/null |sort -S100G -T ./ --parallel=24 -k6,6 -k8,8n > yahs-vs-Csq_v2.0.paf &
```
