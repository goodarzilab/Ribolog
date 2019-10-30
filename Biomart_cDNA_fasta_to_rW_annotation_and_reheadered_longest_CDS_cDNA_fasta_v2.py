#!/usr/bin/python
# scrpit_name.py fasta_in fasta_out annotation_out id_mapper_out <optional> no_x_cds

import sys

##### The difference between this version (v2) and the first code is the order of columns in the ID mapper: transcript ID, gene ID, gene name 

### Step 1: Create a dictionary of genes and the length of their longest CDS

infile = open(sys.argv[1], "r") # fasta of cDNA sequences downloaded from Biomart with header: gene_id | transcript_id | gene_name | cds_start | cds_end | transcript_length
line=infile.readline()

genes_longest_cds={}
i=0

while line:
        if '>' in line:
                i += 1
                lines = line.replace('>','').split('|')
                gene_id = lines[0]
                cds_start = min([ int(x) for x in lines[3].split(';')])
                cds_end = max([ int(x) for x in lines[4].split(';')])
                if len(sys.argv) > 5 and sys.argv[5] == "no_x_cds":
                        cds_end += 3
                l_cds = cds_end - cds_start + 1
                if gene_id not in genes_longest_cds:
                        genes_longest_cds[gene_id] = l_cds
                else:
                        genes_longest_cds[gene_id] = max(l_cds, genes_longest_cds[gene_id])
        line=infile.readline()
    
print("Transcripts examined: " +str(i))
print("Transcripts to be retained: " + str(len(genes_longest_cds)))
#print(genes_longest_cds)
#print(list(genes_longest_cds.items())[:5])
infile.close()
print("Step 1: DONE")

### Step 2: Start over, copy the longest CDS from each gene to the output file.

infile = open(sys.argv[1],"r") # fasta of cDNA sequences downloaded from Biomart with header: gene_id | transcript_id | gene_name | cds_start | cds_end | transcript_length 
outfile1 = open(sys.argv[2], "w") # fasta containing only one transcript per gene, the transcript with the longest CDS
outfile2 = open(sys.argv[3], "w") # Annotation file with columns: transcript, l_tr, l_utr5, l_cds, l_utr3
outfile3 = open(sys.argv[4], "w") # ID mapper file with columns: gene_id, gene_name, transcript 

outfile2.write("transcript\tl_tr\tl_utr5\tl_cds\tl_utr3\n")
outfile3.write("transcript\tgene_id\tgene_name\n")

used_genes = []

line = infile.readline()

while line:
        if '>' in line:
                lines = line.replace('>','').split('|')
                gene_id = lines[0]
                cds_start = min([ int(x) for x in lines[3].split(';')])
                cds_end = max([ int(x) for x in lines[4].split(';')])
                if len(sys.argv) > 5 and sys.argv[5] == "no_x_cds":
                        cds_end += 3
                l_cds = cds_end - cds_start + 1
                if (gene_id not in used_genes) and genes_longest_cds[gene_id] == l_cds:
                        out='y'
                        used_genes.append(gene_id)
                else:
                        out='n'

        if out=='y':
                if '>' in line:
                        tr_id = lines[1]
                        gene_name = lines[2]
                        l_tr = int(lines[-1])
                        l_utr5 = cds_start - 1
                        l_utr3 = l_tr - (l_utr5 + l_cds)
                        line = '>' + tr_id + '\n'
                        line2 = tr_id + '\t' + str(l_tr) + '\t' + str(l_utr5) + '\t' + str(l_cds) + '\t' + str(l_utr3) + '\n'
                        outfile2.write(line2)
                        line3 = tr_id + '\t' + gene_id + '\t' + gene_name +  '\n'
                        outfile3.write(line3)
                
                outfile1.write(line)

        line=infile.readline()

#print(len(used_genes))
infile.close()
outfile1.close()
outfile2.close()
outfile3.close()
print("Step 2: DONE")
