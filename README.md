# Arcobacter WGS data analysis

## Nanopore data

### Basecalling

Three strains seqeunced with Nanopore Minion using the Rapid Barcoding (24) kit and r10.4.1 flow cell.  
Basecalling was done with Dorado v.0.8.2 using all three different models (on Apple M2 Max).  

Fast model as an example: 

```bash 
dorado basecaller \
    fast@latest \
    /Library/MinKNOW/data/Levis_thesis/Arcobacter/20241101_1509_MN45189_FBA32164_2b012041/pod5/ \
    --kit-name SQK-RBK114-24 > Arcobacter_fast.bam

dorado demux \
    --output-dir Arcobacter_fast \
    --no-classify \
    --emit-fastq \
    Arcobacter_fast.bam
```
### Sequence QC and trimming 

Quality control for the data from three different basecalling models to choose the best one.  
_Resources: 1 CPU, 5 Gb mem, 2h_

__fast__ 

```bash 
for barcode in 15 16 17; do
	/projappl/project_2005273/nano_tools/bin/NanoPlot \
  		-o nanoplot_${barcode}_fast \
		-f png \
		--fastq Data/nanopore/Arcobacter_fast/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz
	
	/projappl/project_2005273/nano_tools/bin/nanoQC \
		-o nanoQC_${barcode}_fast \
		Data/nanopore/Arcobacter_fast/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz
done
```

__hac__ 

```bash 
for barcode in 15 16 17; do
	/projappl/project_2005273/nano_tools/bin/NanoPlot \
  		-o nanoplot_${barcode}_hac \
		-f png \
		--fastq Data/nanopore/Arcobacter_hac/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz
	
	/projappl/project_2005273/nano_tools/bin/nanoQC \
		-o nanoQC_${barcode}_hac \
		Data/nanopore/Arcobacter_hac/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz
done
```

__sup__ 

```bash 
for barcode in 15 16 17; do
	/projappl/project_2005273/nano_tools/bin/NanoPlot \
  		-o nanoplot_${barcode}_sup \
		-f png \
		--fastq Data/nanopore/Arcobacter_sup/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz
	
	/projappl/project_2005273/nano_tools/bin/nanoQC \
		-o nanoQC_${barcode}_sup \
		Data/nanopore/Arcobacter_sup/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz
done
```

After evaluating the QC outputs and choosing the best model, reads are trimmed based on the QC results with chopper.  
_Resources: 4 CPU, 5 Gb mem, 30 min_

```bash
best_model="sup"

mkdir trimmed_nanopore

for barcode in 15 16 17; do
	pigz -cd Data/nanopore/Arcobacter_${best_model}/5aead7076fed16715862274e8905bffa4b22d334_SQK-RBK114-24_barcode${barcode}.fastq.gz |\
  	/projappl/project_2005273/nano_tools/bin/chopper --quality 12 --minlength 1000 --headcrop 20 --tailcrop 20 |\
  	pigz > trimmed_nanopore/barcode${barcode}.trimmed.fastq.gz
done
```

### Genome assembly 

The data will be assembled with Flye and Unicycler (miniasm + Racon).  
_Resources: 8 CPU, 20 Gb mem, 4h_

```bash 
for barcode in 15 16 17; do
	/projappl/project_2005273/nano_assembly/bin/flye \
		--nano-raw trimmed_nanopore/barcode${barcode}.fastq.gz \
		--threads $SLURM_CPUS_PER_TASK \
		--out-dir barcode${barcode}_flye

	/projappl/project_2005273/unicyclerbin/unicycler \
		--long trimmed_nanopore/barcode${barcode}.fastq.gz  \
		--threads $SLURM_CPUS_PER_TASK \
		--mode bold \
		--out barcode${barcode}_unicycler
done
```

Assemblies will be compared using Quast.  
_Resources: 1 CPU, 2 Gb mem, 30 min_

```bash
module load quast
quast --output-dir QUAST_nanopore barcode*_flye/assembly.fasta barcode*_unicycler/assembly.fasta 
```

### Genome annotation

When choosing the best assembly methods, resulting genomes are annotated using bakta.  
_Resources: 6 CPU, 10 Gb mem, 2h_

```bash
best_assembler="flye"

for barcode in 15 16 17; do
	/projappl/project_2005273/bakta/bin/bakta \
       barcode${barcode}_${best_assembler}/assembly.fasta \
      --db /scratch/project_2005273/Arcobacter_project/DBs/bakta/db \
      --prefix barcode${barcode} \
      --genus Aliarcobacter \
      --locus barcode${barcode} \
      --threads $SLURM_CPUS_PER_TASK \
      --output barcode${barcode}_bakta
done
```

Make a folder with softlinks for each of the genomes

```bash
mkdir Aliarcobacter_genomes
ln -s barcode*_bakta/*.fna Aliarcobacter_genomes/
```

### Genome completeness estimation  

We will use CheckM2 to estimate the completeness of the assemblies.  
_Resources: 6 CPU, 75 Gb mem, 2h_

```bash
/projappl/project_2005273/tax_tools/bin/checkm2 predict \
      --output-directory CheckM2_out \
      --lowmem \
      --extension .fna \
      --tmpdir $TMPDIR \
      --input Aliarcobacter_genomes/
```

### Taxonomic annotation

We will use GTDB-tk for the taxonomic annotation of the strains.  
_Resources: 6 CPU, 75 Gb mem, 2h_

```bash
export GTDBTK_DATA_PATH="/scratch/project_2005273/Arcobacter_project/DBs/gtdb/release220"

/projappl/project_2005273/tax_tools/bin/gtdbtk classify_wf \
      --out_dir GTDBTK_out \
      --extension .fna \
      --scratch_dir $TMPDIR \
      --tmpdir $TMPDIR \
      --skip_ani_screen \
      --min_perc_aa 0 \
      --pplacer_cpus 1 \
      --cpus $SLURM_CPUS_PER_TASK \
      --genome_dir Aliarcobacter_genomes/
```

