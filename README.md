# Production of fasta sequences from P. ovale whole genome data
<p>An algorithm for producing the whole genome sequence (.fasta format) for <i>Plasmodium ovale</i> species.<br>
 Prior to use this pipeline, we suppose you produced a <i>sorted.bam</i> file at this step (for example using Samtools).<br>
 The algorithm was written for <i><b>Plasmodium ovale</b></i> species, but can also be used for other <i><b>Plasmodium</b></i> species and <b>other organisms.</b></p>
<h3>1. Prerequisites</h3>
 <p>To use this pipeline, you need the following programs:</p>
 <p>- PERL<br>
 - Python 3.x<br>
 - Samtools<br></p>
 <h3>1. Preparing a pileup file for </h3>
 <p> The first step consists in the production of a pileup file using Samtools <i>mpileup</i> that calculates for each position across the genome the depth coverage, the quality and content of the reads:</p>
 <p><code> samtools mpileup -a -f reference_genome.fasta file_sorted.bam > file_sorted.pileup</code></p>
 <p>where reference_genome.fasta is the genome of interest in .fasta format, and file_sorted.bam is the sorted bam file previously produced with samtools <i> sort</i> function.</p>
 <br>
 <h3>2. Generating a table of A, T, G, C and indels content</h3>
 <p>The second step allows to formate a table that indicate the number of A, T, G, C and indels for each position of the genome. For that, we developed a perl script, named extract_info_from_pileup.pl, that takes the previously produced pileup file as an input (option <code>-p</code>).</p>
  <p><code> perl extract_info_from_pileup.pl -p file_sorted.pileup -o ./</code></p>
<p>The algorithm will produce a table, named info_from_pileup.tsv. We then simplify the table to retain only the major bases at each genomic position. Importantly, we used the IUPAC code to produce the major allele, so other letters than A, T, G and C may appear.</p>
 <p><code>cut info_from_pileup.tsv -f 1,2,8 > result.tsv</code></p>







 <p>An example of reference genome in <i>fasta</i> format and corresponding <i>gff</i> file are provided in the <code>data</code> directory.</p>
 <p><code>python3 Scan_gene_coverage.py -p file_coverage.txt -f reference_genome.fasta -g reference_coordinates.gff -o output.txt</code></p>
 <br>
 <h3>3. Citation</h3>
 <p>If you use this program for your own work, please cite:</p>
 <p>Coppée et al. 5WBF: A low-cost and straightforward whole blood filtration method suitable for whole-genome sequencing of <i>Plasmodium falciparum</i> clinical isolates. (2021) In preparation.</i></p>
