Start with human as an example£¬run perl programs£º
1. perl find_all_sgRNA_z_f_c_y.pl hg19_dna.fa out_sgRNA.txt  out_sgRNA_fasta.txt out_sgRNA_gc_t.txt out_nag_fasta.txt out_no_sgRNA.txt
hg19_dna.fa: human genome which can be downloaded from NCBI.
Other output file: out_sgRNA.txt :sgRNA file £» 
out_sgRNA_fasta.txt: sgRNA fasta format; out_nag_fasta.txt: sgRNA fasta format with nag PAM£¨the two files is the input when using bowtie finding the off-target£©¡£
Other files are useless in the sgRNA searching .
2.Running Bowtie or other softwares which uses out_sgRNA_fasta.txt and out_nag_fasta.txt as input file. Off-target within 3 mismatches then can be gotten. The following process is adding the offtarget information into the out_sgRNA.txt The program should be designed by yourself according to your offtarget file and requirement.

