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

Assemblies will be compared using Quast

```bash
module load quast
quast --output-dir QUAST_nanopore barcode*_flye/assembly.fasta barcode*_unicycler/assembly.fasta 
```

### Genome annotation

When choosing the best assembly methods, resulting genomes are annotated using bakta.  

```bash
best_assembler="flye"

for barcode in 15 16 17; do
	/projappl/project_2005273/bakta/bin/bakta \
       barcode${barcode}_${best_assembler}/assembly.fasta \
      --db /projappl/project_2005273/DBs/bakta/db/ \
      --prefix barcode${barcode}  \
      --genus Arcobacter \
      --locus barcode${barcode} \
      --threads $SLURM_CPUS_PER_TASK \
      --output barcode${barcode}_bakta
done
```

Then something else...