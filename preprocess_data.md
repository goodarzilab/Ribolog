# Processing raw Ribosome Profiling sequencing Data

## Requirements:

1. Fastq files for Riboseq data as well as the matched RNA-seq data.
2. A bowtie2 index of the longest_CDS sequences of your organism as well as rRNA/tRNA sequences that can be used to remove contamination. Ribolog provides these files for arabidopsis, fly, human, maize, mouse, rat, worm, yeast and zebrafish.


## Riboseq part

### Trim adapters

```
for f in *.fastq.gz; do
out=${f/fastq.gz/trim.fastq.gz}
echo "cutadapt -j 10 --trimmed-only -a AGATCGGAAGAGCAC -o $out $f"
cutadapt -j 10 -m 15 --trimmed-only -a AGATCGGAAGAGCAC -o $out $f
done
```
### Split Fastq files by multiplexing barcodes

```
zcat riboseq.trim.fastq.gz | /avicenna/hani/anaconda3/bin/fastx_barcode_splitter.pl --bcfile ../barcodes.txt \
--eol --mismatches 1 --prefix 'riboseq' --suffix '.fastq'
```

### Extract UMI barcodes if the protocol includes them:

```
for f in *.c.fastq.gz; do
out=${f/.c.fastq.gz/.c.5cut.trim.fastq.gz}
echo "umi_tools extract --stdin=$f --bc-pattern=NN --log=5processed.log --stdout $out"
umi_tools extract --stdin=$f --bc-pattern=NN --log=5processed.log --stdout $out
done

for f in *.c.5cut.trim.fastq.gz; do
out=${f/.c.5cut.trim.fastq.gz/.c.53cut.trim.fastq.gz}
echo "umi_tools extract --3prime --stdin=$f --bc-pattern=NNNNN --log=3processed.log --stdout $out"
umi_tools extract --3prime --stdin=$f --bc-pattern=NNNNN --log=3processed.log --stdout $out
done

for f in *.c.53cut.trim.fastq.gz; do
out=${f/.c.53cut.trim.fastq.gz/.c.cut2-5.fastq.gz}
echo "zcat $f | perl -n -e'/(@\S+)_(\S\S)_(\S\S\S\S\S)/ ? print ($1."_".$2.$3."\n") : print' | gzip -c > $out";
zcat $f | perl -n -e'/(@\S+)_(\S\S)_(\S\S\S\S\S)/ ? print ($1."_".$2.$3."\n") : print' | gzip -c > $out
done
```

### Remove the rRNA / tRNA contamination reads

```
for f in *.c.cut2-5.fastq.gz; do
out=${f/.c.cut2-5.fastq.gz/.uncontam.fastq.gz}
echo "bowtie2 -p 16 --end-to-end --un-gz=$out -x /avicenna/genomes/mm10/contam/RNAcontam -U $f 2>> stats.txt > aln.out"
bowtie2 -p 16 --end-to-end --un-gz=$out -x /avicenna/genomes/mm10/contam/RNAcontam -U $f 2>> stats.txt > aln.out
done
```

### Align to the longest CDS

```
for f in *.uncontam.fastq.gz; do 
out=${f/.uncontam.fastq.gz/.bam} 
echo "=== $f ===" >> align_stats.txt 
bowtie2 --sensitive --end-to-end -N 1 -p 8 -x references/hg38_longest_CDS -U "$f" 2>> align_stats.txt | samtools sort -@ 4 -m 4G -o "$out" - 
samtools index "$out" 
done
```

## RNA-seq part

### Trim adapters

```
for f in *.fastq.gz; do
out=${f/fastq.gz/trim.fastq.gz}
echo "cutadapt -j 10 --trimmed-only -a AGATCGGAAGAGCAC -o $out $f"
cutadapt -j 10 -m 15 --trimmed-only -a AGATCGGAAGAGCAC -o $out $f
done
```

### Remove the rRNA / tRNA contamination reads

```
for f in *.c.cut2-5.fastq.gz; do
out=${f/.c.cut2-5.fastq.gz/.uncontam.fastq.gz}
echo "bowtie2 -p 16 --end-to-end --un-gz=$out -x /avicenna/genomes/mm10/contam/RNAcontam -U $f 2>> stats.txt > aln.out"
bowtie2 -p 16 --end-to-end --un-gz=$out -x /avicenna/genomes/mm10/contam/RNAcontam -U $f 2>> stats.txt > aln.out
done
```

### Align to the longest CDS

```
for f in *.uncontam.fastq.gz; do 
out=${f/.uncontam.fastq.gz/.bam} 
echo "=== $f ===" >> align_stats.txt 
bowtie2 --sensitive --end-to-end -N 1 -p 8 -x references/hg38_longest_CDS -U "$f" 2>> align_stats.txt | samtools sort -@ 4 -m 4G -o "$out" - 
samtools index "$out" 
done
```