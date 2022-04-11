# Garden Warbler genome annotation readme
Example of genome annotation script with genome, proteins from related-species,
 and without species-specific RNAseq reads

Diagram of what `annotate-vertebrate-no-rnaseq-repeatmodeler.sh` is doing

![annotate-no-rnaseq](https://github.com/jelber2/genomes/blob/main/pics/annotate-no-rnaseq.svg)

`annotate-vertebrate-no-rnaseq-repeatmodeler.sh`

```bash
#! /bin/bash

set -eux

# activate miniconda
. /opt/miniconda3/etc/profile.d/conda.sh

cores2=`echo $cores | awk '{print int(($1/2)+0.5)}'`
cores4=`echo $cores | awk '{print int(($1/4)+0.5)}'`

echo ${cores}
echo ${cores2}
echo ${cores4}
echo ${working_dir}
echo ${prefix}
echo ${genome_path}
echo ${protein_gff}
echo ${pred_gff}
echo ${augustus_species}

cd ${working_dir}

echo -e "\
#-----BLAST and Exonerate Statistics Thresholds\n\
blast_type=ncbi+ #set to 'ncbi+', 'ncbi' or 'wublast'\n\
\n\
pcov_blastn=0.8 #Blastn Percent Coverage Threhold EST-Genome Alignments\n\
pid_blastn=0.85 #Blastn Percent Identity Threshold EST-Genome Aligments\n\
eval_blastn=1e-10 #Blastn eval cutoff\n\
bit_blastn=40 #Blastn bit cutoff\n\
depth_blastn=0 #Blastn depth cutoff (0 to disable cutoff)\n\
\n\
pcov_blastx=0.5 #Blastx Percent Coverage Threhold Protein-Genome Alignments\n\
pid_blastx=0.4 #Blastx Percent Identity Threshold Protein-Genome Aligments\n\
eval_blastx=1e-06 #Blastx eval cutoff\n\
bit_blastx=30 #Blastx bit cutoff\n\
depth_blastx=0 #Blastx depth cutoff (0 to disable cutoff)\n\
\n\
pcov_tblastx=0.8 #tBlastx Percent Coverage Threhold alt-EST-Genome Alignments\n\
pid_tblastx=0.85 #tBlastx Percent Identity Threshold alt-EST-Genome Aligments\n\
eval_tblastx=1e-10 #tBlastx eval cutoff\n\
bit_tblastx=40 #tBlastx bit cutoff\n\
depth_tblastx=0 #tBlastx depth cutoff (0 to disable cutoff)\n\
\n\
pcov_rm_blastx=0.5 #Blastx Percent Coverage Threhold For Transposable Element Masking\n\
pid_rm_blastx=0.4 #Blastx Percent Identity Threshold For Transposbale Element Masking\n\
eval_rm_blastx=1e-06 #Blastx eval cutoff for transposable element masking\n\
bit_rm_blastx=30 #Blastx bit cutoff for transposable element masking\n\
\n\
ep_score_limit=20 #Exonerate protein percent of maximal score threshold\n\
en_score_limit=20 #Exonerate nucleotide percent of maximal score threshold" > maker_bopts.ctl

echo -e "\
#-----Location of Executables Used by MAKER/EVALUATOR\n\
makeblastdb=/usr/bin/makeblastdb #location of NCBI+ makeblastdb executable\n\
blastn=/usr/bin/blastn #location of NCBI+ blastn executable\n\
blastx=/usr/bin/blastx #location of NCBI+ blastx executable\n\
tblastx=/usr/bin/tblastx #location of NCBI+ tblastx executable\n\
formatdb= #location of NCBI formatdb executable\n\
blastall= #location of NCBI blastall executable\n\
xdformat= #location of WUBLAST xdformat executable\n\
blasta= #location of WUBLAST blasta executable\n\
RepeatMasker=/opt/RepeatMasker/RepeatMasker #location of RepeatMasker executable\n\
exonerate=/usr/bin/exonerate #location of exonerate executable\n\
\n\
#-----Ab-initio Gene Prediction Algorithms\n\
snap= #location of snap executable\n\
gmhmme3= #location of eukaryotic genemark executable\n\
gmhmmp= #location of prokaryotic genemark executable\n\
augustus=/genetics/elbers/Augustus/bin/augustus #location of augustus executable\n\
fgenesh= #location of fgenesh executable\n\
tRNAscan-SE= #location of trnascan executable\n\
snoscan= #location of snoscan executable\n\
\n\
#-----Other Algorithms\n\
probuild= #location of probuild executable (required for genemark)" > maker_exe.ctl

echo -e "\
#-----Genome (these are always required)\n\
genome=${genome_path} #genome sequence (fasta file or fasta embeded in GFF3 file)\n\
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic\n\
\n\
#-----Re-annotation Using MAKER Derived GFF3\n\
maker_gff= #MAKER derived GFF3 file\n\
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no\n\
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no\n\
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no\n\
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no\n\
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no\n\
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no\n\
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no\n\
\n\
#-----EST Evidence (for best results provide a file for at least one)\n\
est= #set of ESTs or assembled mRNA-seq in fasta format\n\
altest= #EST/cDNA sequence file in fasta format from an alternate organism\n\
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file\n\
altest_gff= #aligned ESTs from a closly relate species in GFF3 format\n\
\n\
#-----Protein Homology Evidence (for best results provide a file for at least one)\n\
protein= #protein sequence file in fasta format (i.e. from mutiple oransisms)\n\
protein_gff=${protein_gff} #aligned protein homology evidence from an external GFF3 file\n\
\n\
#-----Repeat Masking (leave values blank to skip repeat masking)\n\
model_org= #select a model organism for RepBase masking in RepeatMasker\n\
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker\n\
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner\n\
rm_gff= #pre-identified repeat elements from an external GFF3 file\n\
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no\n\
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)\n\
\n\
#-----Gene Prediction\n\
snaphmm= #SNAP HMM file\n\
gmhmm= #GeneMark HMM file\n\
augustus_species= #Augustus gene prediction species model\n\
fgenesh_par_file= #FGENESH parameter file\n\
pred_gff=${pred_gff} #ab-initio predictions from an external GFF3 file\n\
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)\n\
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no\n\
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no\n\
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no\n\
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs\n\
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no\n\
\n\
#-----Other Annotation Feature Types (features MAKER doesn't recognize)\n\
other_gff= #extra features to pass-through to final MAKER generated GFF3 file\n\
\n\
#-----External Application Behavior Options\n\
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases\n\
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)\n\
\n\
#-----MAKER Behavior Options\n\
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)\n\
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)\n\
\n\
pred_flank=200 #flank for extending evidence clusters sent to gene predictors\n\
pred_stats=1 #report AED and QI statistics for all predictions as well as models\n\
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)\n\
min_protein=30 #require at least this many amino acids in predicted proteins\n\
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no\n\
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no\n\
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no\n\
keep_preds=1 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)\n\
\n\
split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)\n\
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no\n\
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'\n\
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes\n\
\n\
tries=2 #number of times to try a contig if there is a failure for some reason\n\
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no\n\
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no\n\
TMP=/dev/shm #specify a directory other than the system default temporary directory for temporary files" > maker_opts.ctl

echo
echo "Step 1 Get od10_vertebrata"
echo

if test -f "odb10_vertebrata_fasta.tar.gz"
then
  echo "Orthodb must be already downloaded"
  echo
  if test -f "proteins.fasta"
  then
    echo "Orthodb must be already be extracted"
    echo
  else
    tar xvf odb10_vertebrata_fasta.tar.gz
    cat vertebrate/Rawdata/* > proteins.fasta
    rm -fr vertebrate/ odb10_vertebrata_fasta.tar.gz
  fi
else
  wget https://v100.orthodb.org/download/odb10_vertebrata_fasta.tar.gz
  tar xvf odb10_vertebrata_fasta.tar.gz
  cat vertebrate/Rawdata/* > proteins.fasta
  rm -fr vertebrate/ odb10_vertebrata_fasta.tar.gz
fi

echo
echo "Step 2 RepeatModeler"
echo

if test -f "${genome_path%%.*}"-families.fa
then
  echo "RepeatModeler must be already finished"
  echo
  if test -f "${genome_path%%.*}"-families.fa2
  then
    echo "ProtExcluder must be already finished"
    echo
  else
    # Run blastx then ProtExcluder to excluce known protein sequences from RepeatModeler library
    /usr/bin/blastx -num_threads ${cores} -db /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta -evalue 1e-6 \
    -query "${genome_path%%.*}"-families.fa -out "${genome_path%%.*}"-families.fa.blast
    
    /opt/ProtExcluder1.1/ProtExcluder.pl -f 50 "${genome_path%%.*}"-families.fa.blast "${genome_path%%.*}"-families.fa
    mv temp "${genome_path%%.*}"-families.fa2

    rm "${genome_path%%.*}"-families.fa.blast.* "${genome_path%%.*}"-families.fa.ssi "${genome_path%%.*}"-families.fanPr "${genome_path%%.*}"-families.fanoProtFinal
  fi
else
  conda activate repeatmodeler
  conda list
  BuildDatabase -name "${genome_path%%.*}" -engine ncbi "${genome_path%.*}"
  RepeatModeler -pa ${cores4} -engine ncbi -database "${genome_path%%.*}" > repeatmodeler.log 2>&1
  conda deactivate
  if test -f "${genome_path%%.*}"-families.fa2
  then
    echo "ProtExcluder must be already finished"
    echo
  else
    # Run blastx then ProtExcluder to excluce known protein sequences from RepeatModeler library
    /usr/bin/blastx -num_threads ${cores} -db /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta -evalue 1e-6 \
    -query "${genome_path%%.*}"-families.fa -out "${genome_path%%.*}"-families.fa.blast
    
    /opt/ProtExcluder1.1/ProtExcluder.pl -f 50 "${genome_path%%.*}"-families.fa.blast "${genome_path%%.*}"-families.fa
    mv temp "${genome_path%%.*}"-families.fa2

    rm "${genome_path%%.*}"-families.fa.blast.* "${genome_path%%.*}"-families.fa.ssi "${genome_path%%.*}"-families.fanPr "${genome_path%%.*}"-families.fanoProtFinal
  fi
fi


echo
echo "Step 3 RepeatMasker"
echo

if test -f "${genome_path%.*}".masked
then
  echo "RepeatMasker must be already finished"
  cat "${genome_path%.*}".tbl
  echo
else
  conda activate repeatmodeler
  conda list
#  RepeatMasker -pa ${cores4} -e ncbi -species diptera -dir ./ "${genome_path%.*}" -xsmall > repeatmasker.log 2>&1
  RepeatMasker -pa ${cores4} -e ncbi -lib "${genome_path%%.*}"-families.fa2 -dir ./ "${genome_path%.*}" -xsmall > repeatmasker.log 2>&1
  conda deactivate
  cat "${genome_path%.*}".tbl
fi

echo
echo "Step 4 BRAKER"
echo

if test -f ${working_dir}/braker/output.gff3
then
  echo "BRAKER must be already finished"
  cat ${working_dir}/braker/braker.log
  echo
else
  # activate braker2 environment
  conda activate braker2
  conda list

  # set links to programs
  export AUGUSTUS_CONFIG_PATH=/genetics/elbers/Augustus/config/
  export AUGUSTUS_BIN_PATH=/genetics/elbers/Augustus/bin/
  export AUGUSTUS_SCRIPTS_PATH=/genetics/elbers/Augustus/scripts/
  export PROTHINT_PATH=/genetics/elbers/bin/ProtHint-2.5.0/bin/
  export GENEMARK_PATH=/genetics/elbers/bin/gmes_linux_64/
  export GUSHR_PATH=/genetics/elbers/GUSHR/
  export PATH=/genetics/elbers/BRAKER/scripts:$PATH

  # use BRAKER in protein and rnaseq mode
  braker.pl \
  --species=${augustus_species} \
  --prot_seq=proteins.fasta \
  --genome=${genome_path} \
  --augustus_args="--alternatives-from-sampling=true --minexonintronprob=0.2 --minmeanexonintronprob=0.5 --sample=100 --maxtracks=3 --temperature=2" \
  --softmasking --cores=${cores2} > braker2.log.txt 2>&1
  conda deactivate

  cd ${working_dir}/braker

  # after BRAKER2 is finished, convert the hints to MAKER compatible input for proteins
  awk '$3=="CDSpart"' OFS='\t' hintsfile.gff |awk '{$2="protein2genome";$3="match_part";print}' OFS='\t' > maker.proteins.gff3
  paste <(cut -f 1-8 maker.proteins.gff3) \
  <(cut -f 9 maker.proteins.gff3|perl -pe "s/;src=\S+;/;/g"|\
  perl -pe "s/;//g" |perl -pe "s/grp=/ID=/g"|awk '{print $0"_"FNR";";}' OFS='\t') > maker.proteins.gff32

  # now concatenate the augustus and genemark predictions from BRAKER2
  cat augustus.hints.gtf \
  GeneMark-EP/genemark.f.multi_anchored.gtf > braker2.gtf

  # use GFFRead to add missing gene lines to genemark.gtf file, adjust the stop codons,
  # remove transcript with missing start or stop codons
  # sort the output and convert to GFF3 not GTF
  /genetics/elbers/gffread-0.12.4.Linux_x86_64/gffread -g ${genome_path} --adj-stop -J --sort-alpha -E --keep-genes \
  <(perl -pe "s/\t(\S+)(\.t\d+)\n/\ttranscript_id \"\1\2\"; gene_id \"\1\";\n/ if /\ttranscript\t/g" braker2.gtf|perl -pe "s/\t(\S+)\n/\tgene_id \"\1\";\n/ if /\tgene\t/g") -o- 2>gffread.log | \
  perl -pe "s/\ttranscript\t/\tmRNA\t/g" | \
  cat -n - | perl -pe "s/^/ /g" | \
  perl -pe "s/^\s+(\S+)\s+(.+)Parent=(\S+)\n/\2 ID=\3.\1;Parent=\3\n/ if /\tParent/" | \
  perl -pe "s/^\s+\S+\s+(.+)\n/\1\n/ if /^\s+/" > output.gff3
fi

echo
echo "Step 5 MAKER"
echo

if test -f ${working_dir}/braker/maker_gag_output/gag_output/${augustus_species}.genome.proteins.fasta.against.uniprot_trembl_release_2020_06.blast.log
then
  echo "MAKER must be already finished"
  echo "unannotated gene %"
  cat ${working_dir}/braker/unannotated_percent
  cd ${working_dir}/braker/maker_gag_output/gag_output/ && Rscript stats.R && cat genome.stats
  echo
else

  # enable Jean's local PERL library
  cpanm --local-lib=/home/elbersj/perl5 local::lib && eval $(perl -I /home/elbersj/perl5/lib/perl5/ -Mlocal::lib)

  cd ${working_dir}/braker
  # RUN MAKER3 to filter out genes predicted by Augustus and GeneMark that deviate from evidence
  /genetics/elbers/maker-3.01.03/exe/mpich2/bin/mpiexec -n ${cores} /genetics/elbers/maker-3.01.03/bin/maker -base genome ../maker_opts.ctl ../maker_bopts.ctl ../maker_exe.ctl > maker.output.txt 2>&1

  ## first get gff files
  /genetics/elbers/maker-3.01.03/bin/gff3_merge -d genome.maker.output/genome_master_datastore_index.log -n -o genome.all.gff

  ## second get FASTA files
  /genetics/elbers/maker-3.01.03/bin/fasta_merge -d genome.maker.output/genome_master_datastore_index.log -o genome.all.fasta

  ### i Determine the "best" cut-off value for AED (annotation edit distance)
  wget https://raw.githubusercontent.com/mscampbell/Genome_annotation/master/AED_cdf_generator.pl
  perl AED_cdf_generator.pl -b 0.025 genome.all.gff > maker-run1.genome.fasta.aed.cum.frac.below.txt

  ### ii get protein and transcript sequences that have  annotation edit distance < whatever the AED plot in previous step suggests
  cp genome.all.fasta.all.maker.proteins.fasta genome.all.fasta.all.maker.proteins.0-1.0_AED.fasta
  cp genome.all.fasta.all.maker.transcripts.fasta genome.all.fasta.all.maker.transcripts.0-1.0_AED.fasta

  ## third make a id map for all genes in gff file (i.e., instead of maker-gene-01124356 convert to ${prefix}-000000001)
  /genetics/elbers/maker-3.01.03/bin/maker_map_ids --prefix ${prefix} --justify 8 genome.all.gff > genome.all.gff.id.map

  ### i blast high quality maker proteins against uniprot
  /genetics/elbers/diamond-2.0.4.142 \
  blastp --threads ${cores} \
  --ultra-sensitive \
  --max-target-seqs 1 \
  --evalue 1e-6 \
  --outfmt 6 \
  --db /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta \
  --query genome.all.fasta.all.maker.proteins.0-1.0_AED.fasta 1>genome.all.fasta.all.maker.proteins.fasta.all.blast

  ### ii make copies of files (because scripts below overwrite input files)
  awk -F"\t" -v OFS="\t" '$1=="##gff-version 3"||$2=="maker"' genome.all.gff > genome.all.renamed.gff

  cat genome.all.fasta.all.maker.proteins.0-1.0_AED.fasta > genome.all.fasta.all.maker.proteins.renamed.fasta
  cat genome.all.fasta.all.maker.transcripts.0-1.0_AED.fasta > genome.all.fasta.all.maker.transcripts.renamed.fasta
  cat genome.all.fasta.all.maker.proteins.fasta.all.blast > genome.all.fasta.all.maker.proteins.fasta.renamed.blast

  ### iii rename the maker supplied gene names to ids made in the fifth step above
  /genetics/elbers/maker-3.01.03/bin/map_gff_ids genome.all.gff.id.map genome.all.renamed.gff

  ### iv rename the maker supplied protein and transcript names to ids made in the fifth step above
  /genetics/elbers/maker-3.01.03/bin/map_fasta_ids genome.all.gff.id.map genome.all.fasta.all.maker.proteins.renamed.fasta
  /genetics/elbers/maker-3.01.03/bin/map_fasta_ids genome.all.gff.id.map genome.all.fasta.all.maker.transcripts.renamed.fasta

  ### v rename the maker supplied gene names in the BLAST search from sixth step above to ids made in the fifth step above
  /genetics/elbers/maker-3.01.03/bin/map_data_ids genome.all.gff.id.map genome.all.fasta.all.maker.proteins.fasta.renamed.blast

  ### vi annotate the genes using UniProt/Sprot release 2020_05
  /genetics/elbers/maker-3.01.03/bin/maker_functional_gff /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta genome.all.fasta.all.maker.proteins.fasta.renamed.blast genome.all.renamed.gff > genome.all.renamed.annotated.gff

  ### vii retain only gene annotations that are high quality (AED < 0.50)
  cat genome.all.renamed.annotated.gff > genome.all.renamed.annotated.qualityfilter.gff

  ### viii annotate the proteins and transcripts
  /genetics/elbers/maker-3.01.03/bin/maker_functional_fasta /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta genome.all.fasta.all.maker.proteins.fasta.renamed.blast genome.all.fasta.all.maker.proteins.renamed.fasta > genome.all.fasta.all.maker.proteins.renamed.annotated.fasta
  /genetics/elbers/maker-3.01.03/bin/maker_functional_fasta /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta genome.all.fasta.all.maker.proteins.fasta.renamed.blast genome.all.fasta.all.maker.transcripts.renamed.fasta > genome.all.fasta.all.maker.transcripts.renamed.annotated.fasta

  ### ix make a gff file of just the genes (useful for IGV)
  awk -F"\t" -v OFS="\t" '$1=="##gff-version 3"||$3=="gene"' genome.all.renamed.annotated.qualityfilter.gff > genome.all.renamed.annotated.qualityfilter.genes.only.gff

  ### how many genes have no annotations?
  total_genes=`grep -v "^#" genome.all.renamed.annotated.qualityfilter.genes.only.gff|wc -l`
  unannotated_genes=`grep -c "Protein of unknown function" genome.all.renamed.annotated.qualityfilter.genes.only.gff`
  calc $unannotated_genes\/$total_genes\*100 > unannotated_percent

  # should be about <=~10% of total predicted genes

  # convert MAKER annotations into standard format
  python2 /genetics/elbers/genomeannotation-GAG-997e384/gag.py \
  --fasta ${genome_path} \
  --gff genome.all.renamed.annotated.qualityfilter.gff \
  --fix_start_stop \
  --out maker_gag_output > gag_output.log 2>&1

  pigz proteins.fa 
  mkdir -p empty-dir
  rsync --delete -r empty-dir/ genome.maker.output/
  rm -r empty-dir/ genome.maker.output/

  # change directories
  cd maker_gag_output

  pigz genome.fasta

  # Use Annie to add on annotations again in a format recognized by GAG (next step)
  /genetics/elbers/genomeannotation-annie-4bb3980/annie.py \
  -b ../genome.all.fasta.all.maker.proteins.fasta.renamed.blast \
  -g genome.gff \
  -db /genetics/elbers/maker/uniprot_sprot_release_2020_05.fasta > annie.log 2>&1

  # note that if annie finds the same gene it renames both as
  # for example Wing_0 and Wing_1

  # get MAKER gene name and UniProt gene symbol for adding a "Name=" field to the GFF
  cat annie_output.tsv |paste - - |cut -f 1,3 > annie_annotation.txt

  cut -f 1-8 genome.gff |tail -n +2 > genome.col1-8.gff
  cut -f 9 genome.gff |tail -n +2 > genome.col9.gff


  # step to rename in parallel using 75 cores
  rm -fr parallel
  mkdir -p parallel
  cd parallel

  total_lines=`cat ../genome.col9.gff|wc -l`
  num_files=150
  ((lines_per_file = (total_lines + num_files - 1) / num_files))

  split -d --lines=$lines_per_file ../genome.col9.gff

  myfunc(){
  while IFS=$'\t' read a b;do
    perl -pi -e "s/ID=$a\n/ID=$a;Name=$b\n/" x$1
    perl -pi -e "s/Parent=$a\n/Parent=$a;Name=$b\n/" x$1
  done < ../annie_annotation.txt
  }

  export -f myfunc
  seq -w 0 89 > replicates
  seq -w 9000 9059 >> replicates

  # takes about 10-20 minutes
  time parallel --no-notice --jobs ${cores} 'myfunc' :::: replicates > replicates.log 2>&1

  # concatenate the split pieces
  cat x* > ../genome.col9.gene_names.gff

  # go up a directory
  cd ..

  # combine columns 1-8 and 9 to make a new GFF file
  paste genome.col1-8.gff genome.col9.gene_names.gff |sed '1i ##gff-version 3' > genome.reformatted.gff

  # finally, use GAG again to add on the Gene names properly
  # might take up to an hour
  python2 /genetics/elbers/genomeannotation-GAG-997e384/gag.py \
  --fasta ${genome_path} \
  --gff genome.reformatted.gff \
  --anno annie_output.tsv \
  --fix_start_stop \
  --out gag_output > gag_output.log 2>&1

  cd ${working_dir}/braker/maker_gag_output/gag_output

  pigz genome.fasta

  /genetics/elbers/diamond-2.0.4.142 blastp --threads ${cores} --max-target-seqs 1 \
  --db /genetics/elbers/maker/uniprot_trembl_release_2020_06.fasta.gz \
  --query genome.proteins.fasta \
  --outfmt 6 qlen slen \
  --out ${augustus_species}.genome.proteins.fasta.against.uniprot_trembl_release_2020_06.blast \
  > ${augustus_species}.genome.proteins.fasta.against.uniprot_trembl_release_2020_06.blast.log 2>&1

  cp ${augustus_species}.genome.proteins.fasta.against.uniprot_trembl_release_2020_06.blast genome.proteins.fasta.against.uniprot_trembl_release_2020_06.blast

  # final files are genome.* in the directory
  echo ${working_dir}/braker/maker_gag_output/gag_output/

  cat genome.stats
  echo -e 'test <- read.table("genome.proteins.fasta.against.uniprot_trembl_release_2020_06.blast")\ntest2 <- test$V1/test$V2\ntest3 <- test2[test2 >= 0.85]\ntest4 <- test3[test3 <= 1.15]\ntest5 <- test2[test2 < 0.85]\npaste("There are",length(test2),"protein hits, and", length(test4), round(length(test4)/length(test2)*100,2), "%","are between 0.85 and 1.15 (querylength/subjectlength).")\npaste("There are", length(test5), round(length(test5)/length(test2)*100,2),"%", "protein hits less than 0.85 (querylength/subjectlength).")'> stats.R
  Rscript stats.R
  echo "unannotated gene %"
  cat ${working_dir}/braker/unannotated_percent
fi
```

run ``annotate-vertebrate-no-rnaseq-repeatmodeler.sh` with exported variables

```bash
cd /genetics/elbers/warbler2
export cores=75 && \
export working_dir=/genetics/elbers/warbler2/garden-warbler-odb10_vertebrata-no-rnaseq-27-Jan-2020 && \
export prefix=Sybo_ && \
export genome_path=${working_dir}/GCA_014839755.1_bSylBor1.pri_genomic_no_comments_upper.fasta.masked && \
export protein_gff=${working_dir}/braker/maker.proteins.gff32 && \
export pred_gff=${working_dir}/braker/output.gff3 && \
export augustus_species=garden-warbler-odb10_vertebrata-no-rnaseq-27-Jan-2020 && \
bash annotate-vertebrate-no-rnaseq-repeatmodeler.sh | \
> garden-warbler-odb10_vertebrata-no-rnaseq-27-Jan-2020/garden-warbler-odb10_vertebrata-no-rnaseq-27-Jan-2020_annotate-vertebrate-no-rnaseq-repeatmodeler.sh.log 2>&1 &
