# Fly genome assembly readme
```bash
cd /nfs/scistore16/itgrp/jelbers/fly
```

## Step1 - raw HiFi read QC, decontaminate

Note that an alternative to the `sendsketch.sh`
steps is to download
RefSeq bacteria, RefSeq viruses, RefSeq Human
for example then run `seal.sh`
to remove those sequences from the HiFi reads

### Step 1 visual diagram of steps
#### input is ccs.fastq.gz from PacBio CCS program
![decontaminate](https://github.com/jelber2/genomes/blob/main/pics/decontaminate.svg)

`decontaminate-snakefile`

```yaml
#
shell.executable("/bin/bash")

# set number of threads here
THREADS=48
IDS, = glob_wildcards("ccs-reads/{id}.fastq.gz")

# list of rules which are not deployed to slurm
localrules: all, sendsketch, taxids, getLinks, download

rule all:
        input: expand(""raw-reads/{id}-rm-icecream-rm-contam5.fasta.gz", id=IDS)

# use icecreamfinder.sh 38.87
# from BBMap/BBTools https://sourceforge.net/projects/bbmap/ .
# Icecream cones are inverted repeats that
# can appear in PacBio HiFi reads
# icecreamfinder.sh also splits reads at left over
# SMRTBell adapters
rule icecreamfinder:
        input: "ccs-reads/{id}.fastq.gz"
        output: "{id}/decontaminate/{id}-rm-icecream.fasta.gz"
        shell:
            """
            module purge
            module load bbtools/38.87
            icecreamfinder.sh usejni=t ccs=t \
            in={input} \
            out={output} \
            threads={THREADS}
            """

# Use sendsketch.sh 38.87 to assess potential contamination
# from BBMap/BBTools https://sourceforge.net/projects/bbmap/ .
#
# Note that an alternative to the `sendsketch.sh`
# steps is to download
# RefSeq bacteria, RefSeq viruses, RefSeq Human
# for example then run `seal.sh`
# to remove those sequences from the HiFi reads
rule sendsketch:
        input: "{id}/decontaminate/{id}-rm-icecream.fasta.gz"
        output: "{id}/decontaminate/{id}-rm-icecream.fasta.gz.sketch.log"
        shell:
            """
            module purge
            module load bbtools/38.87
            sendsketch.sh \
            in={input} records=1000 > {output}
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
        input: "{id}/decontaminate/genomes/to.download.txt2"
        output: "{id}/decontaminate/genomes.fasta.gz"
        params: "{id}/decontaminate/genomes"
        shell:
            """
            module load bcftools/1.14
            while read i
            do
              wget $i -O - |pigz -dc |bgzip -@4 >> {output}
            done < {input}
            """

### Use seal.sh 38.87
#from BBMap/BBTools https://sourceforge.net/projects/bbmap/ .
#Use seal.sh with minkmerfraction=0.5 (50% 31-mers from reads must match 
#references to be contamination using libraries of "contamination"
#based on sendsketch.sh). seal.sh seems to have better accuracy than kraken2

rule seal:
        input:
            reads="{id}/decontaminate/{id}-rm-icecream.fasta.gz",
            ref="{id}/decontaminate/genomes.fasta.gz"
        output: "raw-reads/{id}-rm-icecream-rm-contam5.fasta.gz"
        shell:
            """
            module purge
            module load bbtools/38.87
            seal.sh threads={THREADS} k=31 ow=t \
            ref=genomes/genomes.fasta.gz \
            minkmerfraction=0.5 \
            in={input} \
            outu={output} \
            out=/dev/null
            """
```

Initialize snakemake conda environment and make a dry run. Extra `perl` code is to make
the labels for the steps easy to interpret
```bash
conda activate snakemake-7.3.6
snakemake --snakefile decontaminate-snakefile --dryrun  --dag all \
> decontaminate.dag
cat decontaminate.dag | \
perl -pe "s/all/HiFi reads ready for assembly/g" | \
perl -pe "s/sendsketch/use sendsketch.sh from BBMap\/BBTools https\:\/\/sourceforge.net\/projects\/bbmap\/ \\\nto make a 'sketch' of the reads \\\nwith some contamination showing up as matches to the sketch/g" | \
perl -pe "s/icecreamfinder\\\nid\: HiFi-reads/remove left over SMRTbell adapters \\\nand remove inverted repeats with icecreamfinder.sh \\\nfrom BBMap\/BBTools https\:\/\/sourceforge.net\/projects\/bbmap\//g" | \
perl -pe "s/seal/use seal.sh from BBMap\/BBTools https\:\/\/sourceforge.net\/projects\/bbmap\/ \\\nto remove reads that have \\\n50% 31-mers matching contamination reference \\\nseal.sh seems to be more accurate than kraken2/g" | \
perl -pe "s/download/use download links to get contamination references\\\nrun dustmasker on each reference/g" | \
perl -pe "s/taxids/get taxon ids from sendsketch output \\\nwith genomes < 35 Mbp as potential bacterial contamination/g" | \
perl -pe "s/getLinks/use NCBI Datasets and taxon ids to get links to contamination references/g" | \
dot -Tsvg > decontaminate.svg
```
actually run snakefile workflow

```bash
snakemake -j 10 --snakefile decontaminate-snakefile \
--printshellcmds --latency-wait 60 --local-cores 1 --cores all --cluster-config cluster.json \
--cluster "sbatch --export=NONE --no-requeue --job-name {cluster.job-name} --mem={cluster.mem} \
--time={cluster.time} --cpus-per-task={cluster.cores} " all > decontaminate-snakefile.log 2>&1 &
```




## Step 2 - assembly of QC'd HiFi reads

### Step 2 visual diagram of steps
Will use peregrine-2021 pg_asm peregrine-r 0.4.13 (main:a95d696) (https://github.com/cschin/peregrine-2021) for assembly 
#### raw here is really the output, decontaminated reads from Step 1
![pg_asm](https://github.com/jelber2/genomes/blob/main/pics/pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads.svg)

`cluster.json`
```json
{
    "__default__":
    {
        "account": "jelbers",
        "mem": "100G",
        "time": "0:30:0",
        "cores" : "48",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_build_sdb1":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_build_sdb2":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_build_sdb3":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_build_sdb4":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_dp_graph4":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_layout4":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "pg_resolve4":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "1",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    }
}
```

Snakefile for assembling the PacBio HiFi reads 
with peregrine-2021 with 1, 2, and 3 rounds of error-correction

`pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads`

```yaml
#
shell.executable("/bin/bash")

# set input path here
IDS, = glob_wildcards("raw-reads/{id}.fasta.gz")

# set number of threads here
THREADS=48
THREADS2=range(1,THREADS+1)
THREADS3=expand(["{threads}"], threads=THREADS2)
THREADS4=[str(item).zfill(2) for item in THREADS3]


# list of rules which are not deployed to slurm
localrules: all, reads1, reads2, reads3, reads4

rule all:
        input:
            input1=expand("{sample}/pg_asm4/layout_m_p.fa", sample=IDS),
            input2=expand("{sample}/pg_asm4/layout_m_a.fa", sample=IDS)
     
# make a file listing the location of the raw HiFi reads
rule reads1:
        input: "raw-reads/{id}.fasta.gz"
        output: "{id}/pg_asm1/reads1"
        shell: """ls `pwd`/{input} > {output}"""

# make a sequence database for the raw HiFi reads
rule pg_build_sdb1:
        input: "{id}/pg_asm1/reads1"
        output:
            database="{id}/pg_asm1/1.seqdb",
            index="{id}/pg_asm1/1.idx"
        params:
            prefix= "{id}/pg_asm1/1"
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                else
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                fi"""

# make a shimmer index for the raw reads
rule pg_build_idx1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx"
#        output: "{id}/pg_asm1/1-001-of-048.dat"
        output: expand(["{id}/pg_asm1/1-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm1/1",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                else
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                fi"""

# overlap the raw reads
rule pg_ovlp1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx",
            dat=expand(["{id}/pg_asm1/1-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm1/1",
            prefix2= "{id}/pg_asm1/overlap1",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""

# error-correct the raw reads
rule pg_ovlp_ec1:
        input:
            sqdb="{id}/pg_asm1/1.seqdb",
            idx="{id}/pg_asm1/1.idx",
            ovlp=expand(["{id}/pg_asm1/overlap1.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm1/overlap1",
            prefix2= "{id}/pg_asm2/2",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""

# make a file listing the locations of the error-corrected reads
rule reads2:
        input: expand(["{id}/pg_asm2/2_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        output: "{id}/pg_asm2/reads2"
        shell: """ls {input} > {output}"""

# make a sequence database for the 1x corrected HiFi reads
rule pg_build_sdb2:
        input: "{id}/pg_asm2/reads2"
        output:
            database="{id}/pg_asm2/2.seqdb",
            index="{id}/pg_asm2/2.idx"
        params:
            prefix= "{id}/pg_asm2/2"
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                else
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                fi"""

# make a shimmer index for the 1x corrected reads
rule pg_build_idx2:
        input:
            sqdb="{id}/pg_asm2/2.seqdb",
            idx="{id}/pg_asm2/2.idx"
        output:expand(["{id}/pg_asm2/2-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm2/2",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                else
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                fi"""

# overlap the 1x corrected reads
rule pg_ovlp2:
        input:
            sqdb="{id}/pg_asm2/2.seqdb",
            idx="{id}/pg_asm2/2.idx",
            dat=expand(["{id}/pg_asm2/2-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm2/overlap2.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm2/2",
            prefix2= "{id}/pg_asm2/overlap2",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""

# error-correct the 1x corrected reads
rule pg_ovlp_ec2:
        input:
            sqdb="{id}/pg_asm2/2.seqdb",
            idx="{id}/pg_asm2/2.idx",
            ovlp=expand(["{id}/pg_asm2/overlap2.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm3/3_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm2/overlap2",
            prefix2= "{id}/pg_asm3/3",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""

# make a file listing the locations of the 2x error-corrected reads
rule reads3:
        input: expand(["{id}/pg_asm3/3_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        output: "{id}/pg_asm3/reads3"
        shell: """ls {input} > {output}"""
        
# make a sequence database for the 2x corrected HiFi reads
rule pg_build_sdb3:
        input: "{id}/pg_asm3/reads3"
        output:
            database="{id}/pg_asm3/3.seqdb",
            index="{id}/pg_asm3/3.idx"
        params:
            prefix= "{id}/pg_asm3/3"
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                else
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                fi"""

# make a shimmer index for the 2x corrected reads
rule pg_build_idx3:
        input:
            sqdb="{id}/pg_asm3/3.seqdb",
            idx="{id}/pg_asm3/3.idx"
        output: expand(["{id}/pg_asm3/3-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm3/3",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                else
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                fi"""

# overlap the 2x corrected reads
rule pg_ovlp3:
        input:
            sqdb="{id}/pg_asm3/3.seqdb",
            idx="{id}/pg_asm3/3.idx",
            dat=expand(["{id}/pg_asm3/3-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm3/overlap3.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm3/3",
            prefix2= "{id}/pg_asm3/overlap3",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""

# error-correct the 2x corrected reads
rule pg_ovlp_ec3:
        input:
            sqdb="{id}/pg_asm3/3.seqdb",
            idx="{id}/pg_asm3/3.idx",
            ovlp=expand(["{id}/pg_asm3/overlap3.{threads2}"], threads2=THREADS4, allow_missing=True)
        output: expand(["{id}/pg_asm4/4_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm3/overlap3",
            prefix2= "{id}/pg_asm4/4",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp_ec {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""

# make a file listing the locations of the 3x error-corrected reads
rule reads4:
        input: expand(["{id}/pg_asm4/4_{threads2}.fa"], threads2=THREADS4, allow_missing=True)
        output: "{id}/pg_asm4/reads4"
        shell: """ls {input} > {output}"""
        
# make a sequence database for the 3x corrected HiFi reads
rule pg_build_sdb4:
        input: "{id}/pg_asm4/reads4"
        output:
            database="{id}/pg_asm4/4.seqdb",
            index="{id}/pg_asm4/4.idx"
        params:
            prefix= "{id}/pg_asm4/4"
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                else
                    pg_build_sdb --input {input} --out_prefix {params.prefix}
                fi"""

# make a shimmer index for the 3x corrected reads
rule pg_build_idx4:
        input:
            sqdb="{id}/pg_asm4/4.seqdb",
            idx="{id}/pg_asm4/4.idx"
        output: expand(["{id}/pg_asm4/4-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm4/4",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                else
                    pg_build_idx -k {params.k} -r {params.r} -w {params.w} {input.sqdb} {input.idx} {params.prefix} {params.chunks} {params.chunks}
                fi"""

# overlap the 3x corrected reads
rule pg_ovlp4:
        input:
            sqdb="{id}/pg_asm4/4.seqdb",
            idx="{id}/pg_asm4/4.idx",
            dat=expand(["{id}/pg_asm4/4-0{threads2}-of-0{threads}.dat"], threads2=THREADS4, threads=THREADS, allow_missing=True)
        output: expand(["{id}/pg_asm4/overlap4.{threads2}"], threads2=THREADS4, allow_missing=True)
        params:
            prefix1= "{id}/pg_asm4/4",
            prefix2= "{id}/pg_asm4/overlap4",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                else
                    pg_ovlp {input.sqdb} {input.idx} {params.prefix1} {params.prefix2} {params.chunks} {params.chunks}
                fi"""
        
# make a sequence graph of the 3x-corrected reads
rule pg_dp_graph4:
        input: expand(["{id}/pg_asm4/overlap4.{threads2}"], threads2=THREADS4, allow_missing=True)
        output:
            multiext("{id}/pg_asm4/overlap4", "_g.dat", "_utg.gfa", "_utg.dat", "_utg0.dat", "_g0.dat", "_layout.dat")
        params:
            b="6",
            prefix1= "{id}/pg_asm4/overlap4"
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_dp_graph -b {params.b} --prefix {params.prefix1} --out_prefix {params.prefix1}
                else
                    pg_dp_graph -b {params.b} --prefix {params.prefix1} --out_prefix {params.prefix1}
                fi"""

# from the 3x-corrected-read-graph-layout, generate contigs
rule pg_layout4:
        input:
            sqdb="{id}/pg_asm4/4.seqdb",
            idx="{id}/pg_asm4/4.idx",
            layout="{id}/pg_asm4/overlap4_layout.dat"
        output:
            primary="{id}/pg_asm4/layout_m.fa",
            alternate="{id}/pg_asm4/layout_e0.fa"
        params:
            prefix1= "{id}/pg_asm4/layout",
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_layout {input.sqdb} {input.idx} --layout_file {input.layout} --output {params.prefix1}
                else
                    pg_layout {input.sqdb} {input.idx} --layout_file {input.layout} --output {params.prefix1}
                fi"""

# from 3x-corrected reads generated contigs resolve haplotypes for the contigs to primary
# layout_m_p.fa and alternate layout_m_a.fa
rule pg_resolve4:
        input: "{id}/pg_asm4/layout_m.fa"
        output:
            multiext("{id}/pg_asm4/layout_m", "_p.fa", "_a.fa")
        params:
            k="56",
            r="6",
            w="80",
            prefix= "{id}/pg_asm4/layout_m",
            chunks= THREADS
        shell: """module purge
                if [ `which pg_asm |wc -l` -eq 0 ]
                then
                    export "PATH=/nfs/scistore16/itgrp/jelbers/git/peregrine-2021/target/release:$PATH"
                    pg_resolve -k {params.k} -r {params.r} -w {params.w} -f {input} -o {params.prefix}
                else
                    pg_resolve -k {params.k} -r {params.r} -w {params.w} -f {input} -o {params.prefix}
                fi"""
```

Initialize snakemake conda environment and make a dry run. Extra `perl` code is to make
the labels for the steps easy to interpret
```bash
conda activate snakemake-7.3.6
snakemake --snakefile pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads --dryrun  --dag all \
> pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads.dag
cat pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads.dag | \
perl -pe "s/reads1\\\nid\: HiFi-reads-rm-icecream-rm-contam5/make a file listing the location of the raw HiFi reads/g" | \
perl -pe "s/pg_build_sdb1/make a sequence database for the raw HiFi reads/g"| \
perl -pe "s/pg_build_idx1/make a shimmer index for the raw reads/g"| \
perl -pe "s/pg_ovlp1/overlap the raw reads/g"| \
perl -pe "s/pg_ovlp_ec1/error-correct the raw reads/g"| \
perl -pe "s/reads2/make a file listing the locations of the error-corrected reads/g"| \
perl -pe "s/pg_build_sdb2/make a sequence database for the 1x corrected HiFi reads/g"| \
perl -pe "s/pg_build_idx2/make a shimmer index for the 1x corrected reads/g"| \
perl -pe "s/pg_ovlp2/overlap the 1x corrected reads/g"| \
perl -pe "s/pg_ovlp_ec2/error-correct the 1x corrected reads/g"| \
perl -pe "s/reads3/make a file listing the locations of the 2x error-corrected reads/g"| \
perl -pe "s/pg_build_sdb3/make a sequence database for the 2x corrected HiFi reads/g"| \
perl -pe "s/pg_build_idx3/make a shimmer index for the 2x corrected reads/g"| \
perl -pe "s/pg_ovlp3/overlap the 2x corrected reads/g"| \
perl -pe "s/pg_ovlp_ec3/error-correct the 2x corrected reads/g"| \
perl -pe "s/reads4/make a file listing the locations of the 3x error-corrected reads/g"| \
perl -pe "s/pg_ovlp4/overlap the 3x corrected reads/g"| \
perl -pe "s/pg_dp_graph4/make a sequence graph of the 3x-corrected reads/g"| \
perl -pe "s/pg_layout4/from the 3x-corrected-read-graph-layout, generate contigs/g"| \
perl -pe "s/pg_resolve4/from 3x-corrected reads generated contigs make primary and alternate contigs/g" | \
perl -pe "s/pg_build_sdb4/make a sequence database for the 3x corrected HiFi reads/g"| \
perl -pe "s/pg_build_idx4/make a shimmer index for the 3x corrected reads/g" | \
perl -pe "s/all/ready for purge_dups/g" | \
dot -Tsvg > pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads.svg
```

actually run snakefile workflow

```bash
snakemake -j 10 --snakefile pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads \
--printshellcmds --latency-wait 60 --local-cores 1 --cores all --cluster-config cluster.json \
--cluster "sbatch --export=NONE --no-requeue --job-name {cluster.job-name} --mem={cluster.mem} \
--time={cluster.time} --cpus-per-task={cluster.cores} " all > pg_asm-1_2_3_rounds-ec-resolve_snakefile-threads.log 2>&1 &
```

## Step 3 run purge_dups

### Step 3 visual diagram of steps
![purge_dups](https://github.com/jelber2/genomes/blob/main/pics/purge_dups.svg)

`purge_dups-snakefile`

```yaml
#
shell.executable("/bin/bash")

# set number of threads here
IDS, = glob_wildcards("raw-reads/{id}.fasta.gz")

THREADS=48
THREADS2=range(1,THREADS+1)
THREADS3=expand(["{threads}"], threads=THREADS2)
THREADS4=[str(item).zfill(2) for item in THREADS3]
READS=expand("{id}/pg_asm4/4_{threads2}.fa", threads2=THREADS4, allow_missing=True)
REFERENCE="{id}/pg_asm4/layout_m_p.fa"

# list of rules which are not deployed to slurm
localrules: all, pbcstat, calcuts, split_fa, purge_dups, get_seqs

rule all:
        input:
            purged=expand("{id}/pg_asm4/test.purged.fa", id=IDS),
            hap=expand("{id}/pg_asm4/test.hap.fa", id=IDS)

# map reads to assembly to purge with minimap2
rule minimap2_reads:
        input:
            reads=READS,
            ref=REFERENCE
        output: "{id}/pg_asm4/4.paf.gz"
        shell: """
            module purge
            module load compression-tools/20220329
            module load minimap2/2.24
            minimap2 -t {THREADS} --secondary=no -k19 -w10 -O5,56 -E4,1 -A2 -B5 -z400,50 -r2000 \
            -g5000 {input.ref} {input.reads} | pigz -p {THREADS} -c - > {output}
            """

# get coverage
rule pbcstat:
        input: "{id}/pg_asm4/4.paf.gz"
        output:
            multiext("{id}/pg_asm4/PB.","stat", "base.cov", "cov.wig")
        params: "{id}/pg_asm4"
        shell:
            """
            module purge
            if [ `which purge_dups |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/purge_dups/src:$PATH"
                pbcstat {input} -O {params}
            else
                pbcstat {input} -O {params}
             fi
             """

# get cutoffs
rule calcuts:
        input: "{id}/pg_asm4/PB.stat"
        output: "{id}/pg_asm4/PB.cutoffs"
        shell:
            """
            module purge
            if [ `which purge_dups |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/purge_dups/src:$PATH"
                calcuts {input} > {output}
            else
                calcuts {input} > {output}        
             fi
             """

# split reference
rule split_fa:
        input: REFERENCE
        output: "{id}/pg_asm4/split.fa"
        shell:
            """
            module purge
            if [ `which purge_dups |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/purge_dups/src:$PATH"
                split_fa {input} > {output}
            else
                split_fa {input} > {output}
             fi
             """

# map split assembly to itself with minimap2
rule minimap2_ref:
        input:
            ref="{id}/pg_asm4/split.fa"
        output: "{id}/pg_asm4/split.fa.paf.gz"
        shell: """
            module purge
            module load compression-tools/20220329
            module load minimap2/2.24
            minimap2 -t {THREADS} -xasm5 -DP \
            {input.ref} {input.ref} | pigz -p {THREADS} -c - > {output}
            """

# purge_dups
rule purge_dups:
        input:
            coverage="{id}/pg_asm4/PB.base.cov",
            cutoffs="{id}/pg_asm4/PB.cutoffs",
            split_paf="{id}/pg_asm4/split.fa.paf.gz"
        output: "{id}/pg_asm4/PB.bed"
        shell:
            """
            module purge
            if [ `which purge_dups |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/purge_dups/src:$PATH"
                purge_dups -2 -T {input.cutoffs} -c {input.coverage} {input.split_paf} > {output}
            else
                purge_dups -2 -T {input.cutoffs} -c {input.coverage} {input.split_paf} > {output}
             fi
             """

# actually do purging
rule get_seqs:
        input:
            bed="{id}/pg_asm4/PB.bed",
            ref=REFERENCE
        output:
            purged="{id}/pg_asm4/test.purged.fa",
            hap="{id}/pg_asm4/test.hap.fa"
        params: "{id}/pg_asm4/test"
        shell:
            """
            module purge
            if [ `which purge_dups |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/purge_dups/src:$PATH"
                get_seqs -p {params} {input.bed} {input.ref}
            else
                get_seqs -p {params} {input.bed} {input.ref}  
             fi
             """
```

Initialize snakemake conda environment and make a dry run. Extra `perl` code is to make
the labels for the steps easy to interpret
```bash
conda activate snakemake-7.3.6
snakemake --snakefile purge_dups-snakefile --dryrun  --dag all \
> purge_dups.dag
cat purge_dups.dag | \
perl -pe "s/minimap2_reads\\\nid\: HiFi-reads-rm-icecream-rm-contam5/map 3x error-corrected reads to primary contigs/g" | \
perl -pe "s/pbcstat/compute read coverage/g"| \
perl -pe "s/calcuts/compute coverage cut-offs/g"| \
perl -pe "s/purge_dups/predict duplicate regions/g"| \
perl -pe "s/get_seqs/do actual purging/g"| \
perl -pe "s/split_fa\\\nid\: HiFi-reads-rm-icecream-rm-contam5/split primary contigs/g" | \
perl -pe "s/minimap2_ref/map split primary contigs to itself/g" | \
perl -pe "s/all/ready for k-mer evaluation/g" | \
dot -Tsvg > purge_dups.svg
```

actually run snakefile workflow

```bash
snakemake -j 10 --snakefile purge_dups-snakefile \
--printshellcmds --latency-wait 60 --local-cores 4 --cores all \
--cluster "sbatch --export=NONE --no-requeue --job-name {rule} --mem=64g \
--time=1:00:00 --cpus-per-task=48 " all > purge_dups-snakefile.log 2>&1 &
```

## Step 4 evaluate purged assembly with k-mers

run FastK (https://github.com/thegenemyers/FASTK) 
and MerquryFK (https://github.com/thegenemyers/MERQURY.FK)

### Step 4 visual diagram of steps
![merquryfk](https://github.com/jelber2/genomes/blob/main/pics/merquryfk.svg)


`merquryfk-cluster.json`

```json
{
    "__default__":
    {
        "account": "jelbers",
        "mem": "100G",
        "time": "0:30:0",
        "cores" : "4",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    },
    "fastk":
    {
        "account": "jelbers",
        "time": "0:30:0",
        "mem": "100G",
        "cores" : "48",
        "job-name" : "{rule}",
        "partition" : "defaultp"
    }
}
```

`merquryfk-snakefile`

```yaml
#
shell.executable("/bin/bash")

# set reference directory
IDS, = glob_wildcards("raw-reads/{id}.fasta.gz")

# set number of threads here
THREADS=48
THREADS2=range(1,THREADS+1)
THREADS3=expand(["{threads}"], threads=THREADS2)
THREADS4=[str(item).zfill(2) for item in THREADS3]
READS1=expand(["{id}/pg_asm4/4_{threads2}.fa"], id=IDS, threads2=THREADS4)
REFERENCE1=expand(["{id}/pg_asm4/test.purged.fa"], id=IDS)

READS2=expand(["{id}/pg_asm4/4_{threads2}"], id=IDS, threads2=THREADS4)
REFERENCE2="test.purged"

# list of rules which are not deployed to slurm
localrules: all

rule all:
        input:
            input4=expand("{sample}/pg_asm4/merqury.fk.completeness.stats", sample=IDS)
     
# run fastk
rule fastk:
        input: READS1
        output:
            multiext("{id}/pg_asm4/reads.", "prof", "ktab", "hist")
        params: "{id}/pg_asm4/"
        shell:
            """
            module purge
            if [ `which FastK |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/FASTK-1.0:$PATH"
                FastK -P{params} -T{THREADS} -t1 -k56 -p -N{params}reads -v {READS2}
            else
                FastK -P{params} -T{THREADS} -t1 -k56 -p -N{params}reads -v {READS2}
             fi
             """

# run merqury.fk
rule merquryfk:
        input: 
            multiext("{id}/pg_asm4/reads.", "prof", "ktab", "hist"),
            REFERENCE1
        output:
            multiext("{id}/pg_asm4/merqury.fk.", "qv", "spectra-asm.fl.pdf", "completeness.stats")
        params: "{id}/pg_asm4/"
        shell:
            """
            module purge
            module load R/4.1.2
            if [ `which MerquryFk |wc -l` -eq 0 ]
            then
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/MERQURY.FK:$PATH"
                export "PATH=/nfs/scistore16/itgrp/jelbers/bin/FASTK-1.0:$PATH"
                cd {params}
                MerquryFK -v -f -pdf -P./ reads {REFERENCE2} merqury.fk
            else
                cd {params}
                MerquryFK -v -f -pdf -P./ reads {REFERENCE2} merqury.fk
             fi
             """
```

Initialize snakemake conda environment and make a dry run. Extra `perl` code is to make
the labels for the steps easy to interpret
```bash
conda activate snakemake-7.3.6
snakemake --snakefile merquryfk-snakefile --dryrun  --dag all \
> merquryfk.dag
cat merquryfk.dag | \
perl -pe "s/fastk\\\nid\: HiFi-reads-rm-icecream-rm-contam5/count 56-mers from 3x error-corrected reads with FastK/g" | \
perl -pe "s/merquryfk/run MerquryFK with counted 56-mers and purged primary contigs/g"| \
perl -pe "s/all/ready to evaluate quality value of assembly,etc/g" | \
dot -Tsvg > merquryfk.svg
```



actually run snakefile workflow

```bash
snakemake -j 10 --snakefile merquryfk-snakefile \
--printshellcmds --latency-wait 60 --local-cores 1 --cores all --cluster-config merquryfk-cluster.json \
--cluster "sbatch --export=NONE --no-requeue --job-name {cluster.job-name} --mem={cluster.mem} \
--time={cluster.time} --cpus-per-task={cluster.cores} " all > merquryfk-snakefile.log 2>&1 &
```





MerquryFK suggests the following

```bash
cat HiFi-reads-rm-icecream-rm-contam5/pg_asm4/merqury.fk.completeness.stats | \
column -ts $'\t'

Assembly     Region  Found      Total      % Covered
test.purged  all     998492792  998510479  100.00

# note that Found and Total are the number of 56-mers

cat HiFi-reads-rm-icecream-rm-contam5/pg_asm4/merqury.fk.qv | \
column -ts $'\t'

Assembly     No Support  Total      Error %  QV
test.purged  17687       762435740  0.0000   63.8 # estimated quality value

# No Support and Total are values in bases
```

![spectra-plot](https://github.com/jelber2/genomes/blob/main/pics/merqury.fk.test.purged.spectra-cn.fl.svg)

Copy number spectra plot of 56-mers against the genome
showing very high heterozygosity and good collapsing of haplotypes

# BUT HOW ACCURATE ARE MERQURY.FK RESULTS WITHOUT ADDITIONAL DATA? SEE [simulation.md](https://github.com/jelber2/genomes/blob/main/simulation.md)
