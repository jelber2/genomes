# Analysis of Chrysomallon squamiferum

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
snakemake -j 10 --snakefile decontaminate-snakefile --printshellcmds --latency-wait 60 --local-cores 4 --cores all --cluster "sbatch --export=NONE --no-requeue --job-name {rule} --mem=100g \
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
~/bin/rasusa-0.6.1-x86_64-unknown-linux-musl/rasusa --seed 1 --coverage 50 --genome-size 400M --input SRR12763791_1.fastq.gz \
2> rasusa.log |seqtk seq -A > SRR12763791.50x.fasta
```

## Quality and adapter trimming
bbduk 38.82

bbduk.sh threads=24 in1=SRR8599727_decon_1.fastq.gz in2=SRR8599727_decon_2.fastq.gz \
out1=SRR8599727_decon_trim_1.fastq.gz out2=SRR8599727_decon_trim_2.fastq.gz \
ref=~/bin/bbmap-38.94/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=15 > bbduk2.log 2>&1 &

## Error-correct decontaminated, QC'd Illumina reads
tadpole.sh in1=SRR8599727_decon_trim_1.fastq.gz in2=SRR8599727_decon_trim_2.fastq.gz \
out1=<corrected> mode=correct


## Use GraphAligner commit # cf7f5db and decontaminated, QC'd, and error-corrected Illumina WGS reads to correct 50x Nanopore
/nfs/scistore16/itgrp/jelbers/bin/bcalm/build/bcalm -nb-cores 40 -in filenames -kmer-size 67 -abundance-min 3 > tmp/bcalm_stdout.txt 2> tmp/bcalm_stderr.txt
  
  
```sh
git show

commit cf7f5dbaab1e852d4a2e0c5834f56c1b49424b04 (HEAD -> master, tag: v1.0.16-osx, origin/master, origin/HEAD)
Author: Mikko Rautiainen <m_rautiainen@hotmail.com>
Date:   Fri Mar 25 15:03:28 2022 -0400

    correct mxm library in readme
```

# https://github.com/maickrau/GraphAligner
GraphAlignerPath: /nfs/scistore16/itgrp/jelbers/bin/GraphAligner/bin/GraphAligner
# https://github.com/GATB/bcalm
BcalmPath: /nfs/scistore16/itgrp/jelbers/bin/bcalm/build/bcalm
# https://github.com/GATB/bcalm/blob/master/scripts/convertToGFA.py
BcalmConvertPath: /nfs/scistore16/itgrp/jelbers/bin/bcalm/scripts/convertToGFA.py
# https://github.com/mourisl/Lighter

  
  ### make lower case bases upper case and put on a single line

```sh
seqtk seq -Ul0 output/corrected.50x.fa > SRR12763791.50x.GraphAligner.fasta
```


