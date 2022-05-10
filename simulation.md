# Simulation based on human chromosome 8 and decontaminated reads from ![fly](https://github.com/jelber2/genomes/blob/main/fly-genome.md)
```bash
cd /nfs/scistore16/itgrp/jelbers/chr8
```

## Step1 Simulate Chr8 reads with ![HI.SIM](https://github.com/thegenemyers/HI.SIM)

### Get Telomere to Teleomere human chr8 assembly
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/894/425/GCA_016894425.1_ASM1689442v1/GCA_016894425.1_ASM1689442v1_genomic.fna.gz
# convert lowercase to uppercase bases
~/bin/seqtk/seqtk seq -U GCA_016894425.1_ASM1689442v1_genomic.fna.gz > GCA_016894425.1_ASM1689442v1_genomic_upper.fasta
```

### Put program binaries in PATH
```bash
export PATH="/nfs/scistore16/itgrp/jelbers/bin/FASTK-1.0:$PATH"
export PATH="/nfs/scistore16/itgrp/jelbers/bin/MERQURY.FK:$PATH"
export PATH="/nfs/scistore16/itgrp/jelbers/bin/HI.SIM:$PATH"
```

### Count 40-mers with FastK and make them symmetric with Symmex
```bash
THREADS=20
FastK -k40 -T${THREADS} -P./ -N./Hifi -t1 -p -v HiFi-reads-rm-icecream-rm-contam5.fasta.gz > FastK.log 2>&1 &
Symmex -P./ -v -T${THREADS} Hifi Hifi2 > Symmex.log 2>&1 &

### this step is to rename files produced by FastK to those needed by HI.SIM
```bash 
for i in `seq 1 ${THREADS}`
do
  cp .Hifi.pidx.${i} .Hifi2.pidx.${i}
  cp .Hifi.prof.${i} .Hifi2.prof.${i}
done

cp Hifi.prof Hifi2.prof
cp Hifi.hist Hifi2.hist
```

### Run HImodel and HIsim
```bash
HImodel -v -T${THREADS} -v -g2:51 -e2 -oHifi2 Hifi2 > HImodel.log 2>&1 &
HIsim GCA_016894425.1_ASM1689442v1_genomic_upper Hifi2 -oreference -h -e \
-f -c50 -m10000 -s1000 -x1000 -v -C -U -r1 > HIsim.log 2>&1 &
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
cat reference/pg_asm4/merqury.fk.completeness.stats | \
column -ts $'\t'

Assembly     Region  Found      Total      % Covered
test.purged  all     139289577  139289682  100.00

# note that Found and Total are the number of 56-mers

cat reference/pg_asm4/merqury.fk.qv | \
column -ts $'\t'

Assembly     No Support  Total      Error %  QV
test.purged  105         146068892  0.0000   78.9 # estimated quality value

# No Support and Total are values in bases
```

# Ok, so that is using the reads used for assembly to assess quality.
# What if we use perfect reads?

use BBTools/BBMap 38.82 (https://sourceforge.net/projects/bbmap/)
to simulate 9000-12000 bp reads at 15x coverage that have no errors against the original chr8 assembly
```bash
randomreads.sh ow=t seed=1 ref=GCA_016894425.1_ASM1689442v1_genomic_upper.fasta illuminanames=t \
pacbio=t pbmin=0 pbmax=0 coverage=15 paired=f gaussianlength=t minlength=9000 midlength=10000 \
maxlength=12000 out=merquryfk/perfect-hifi.fasta.gz 2>/dev/null &
```

Use FastK to count 56-mers
```bash
cd merquryfk
FastK -P./ -T34 -t1 -k56 -p -Nreads -v perfect-hifi
```

Run MerquryFK to assess quality of de novo assembly using perfect reads
```bash
ln -s ../reference/pg_asm4/test.purged.fa pg_asm3x-purge_dups.fasta
MerquryFK -v -f -pdf -P./ reads pg_asm3x-purge_dups pg_asm3x-purge_dups
```

```bash
# completeness
column -ts $'\t' pg_asm3x-purge_dups.completeness.stats

Assembly             Region  Found      Total      % Covered
pg_asm3x-purge_dups  all     139282167  139282896  100.00

# quality value
column -ts $'\t' pg_asm3x-purge_dups.qv

Assembly             No Support  Total      Error %  QV
pg_asm3x-purge_dups  731         146068892  0.0000   70.5
```

Instead of MerquryFK and FastK, let's try [YAK](https://github.com/lh3/yak)

```bash
# count 31-mers and discard singletons
~/bin/yak/yak count -b37 -t32 -o perfect-hifi.fasta.yak perfect-hifi.fasta 2>/dev/null &

# quality value estimation
yak qv -t34 perfect-hifi.fasta.yak pg_asm3x-purge_dups.fasta > pg_asm3x-purge_dups.fasta.yak.qv.txt 2>/dev/null &

# quality values and explanation
head -n 4 pg_asm3x-purge_dups.fasta.yak.qv.txt && tail -n 4 pg_asm3x-purge_dups.fasta.yak.qv.txt

CC	CT  kmer_occurrence    short_read_kmer_count  raw_input_kmer_count  adjusted_input_kmer_count
CC	FR  fpr_lower_bound    fpr_upper_bound
CC	ER  total_input_kmers  adjusted_error_kmers
CC	CV  coverage
FR	0	1
ER	146068942	2737.287
CV	1.000
QV	58.510	62.186

So QV is predicted between 58.510-62.186
completeness/coverage is predicted to be 100 %
```
