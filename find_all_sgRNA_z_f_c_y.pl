#!/usr/bin/perl -w
#find all possible sgRNA  binding sites.

my ($infile1,$outfile1,$outfile2,$outfile3,$outfile4,$outfile5)= @ARGV;
open (INA, $infile1) || die;      #mm10.fa
open (OUTA, ">$outfile1") || die;   #sgRNA 20bp
open (OUTB, ">$outfile2") || die;  #sgRNA fastq
open (OUTC, ">$outfile3") || die;  #gc and poly t
open (OUTD, ">$outfile4") || die;  #sgRNA nag fastq
open (OUTE, ">$outfile5") || die;  #sgRNA number of each chromatin


my $seq = "";
my $no_sgRNA=0;
my $bp23="";
my $bp20="";
my $bp21="";
my $polyt=0;
my $escore=0;
my @chromosome;
my $id_chr=1;
my @no_sgRNA_chr;
#my $chr_flag=0;
my $bp23_ag="";
$no_sgRNA_chr[0]=0;
while (<INA>) {
	chomp;
	
	if($_ =~ />chr(\w+)$/) {
		#$chr_flag++;
		$chromosome[$id_chr]=$1;
		$id_chr++;
		if((length($seq) !=0)&&($id_chr>2)) {
			my $no_gg=0;	
			
			my $flag=0;
			my $start_bp=21;
			# find sgRNA
			while($flag !=(-1)) {
				$flag = index($seq, "GG",$start_bp);
				my $start_23=0;
				if($flag !=(-1))  {$start_23=$flag-21;}
				$location=$start_23;
				$bp23 = substr ($seq,$start_23,23);
				$start_bp=$flag+1;
				
					$no_sgRNA ++;
					$bp20 = substr($bp23,0,20); 
					$bp21 = substr($bp23,0,21); 
					$bp23_ag=$bp21 ."AG";
					
					my $count1=$bp20=~tr/G/G/;
					my $count2=$bp20=~tr/C/C/;
					my $gc = ($count1+$count2)/20;
					$escore=0;
					if(($gc<0.75)&&($gc>0.25)) {
						$escore= $escore+5;
					}
					if($bp20 =~ /TTTT/) {
						;
					}
					else {
						$escore= $escore+5;
					}
					
					printf OUTB ">$no_sgRNA\n$bp23\n";
					printf OUTD ">$no_sgRNA\n$bp23_ag\n";
					
					printf OUTA "$bp20\t$chromosome[$id_chr-2]\t$location\t-\t$escore\n"  ;
					
			}  ## while
			 $start_bp=0;
			$flag=0;
			$seq =$seq."MMMMMMMMMMMMMMMMMMMMMM";
			while($flag !=(-1)) {
				$flag = index($seq, "CC",$start_bp);
			
				my $start_23=0;
				if($flag !=(-1))  {$start_23=$flag;}
				$location=$start_23+3;
				$bp23 = substr ($seq,$start_23,23);
				
				$start_bp=$flag+1;
				
					#$no_sgRNA ++;
				if(substr($bp23,22,1) ne "M") {	
					$no_sgRNA ++;
				#	printf "a\n";
					$bp20 = substr($bp23,3,20); 
					
					$bp21 = substr($bp23,2,21); 
					
					$bp20= reverse $bp20;
					$bp20=~tr/ACGT/TGCA/;
					$bp21= reverse $bp21;
					$bp21=~tr/ACGT/TGCA/;
					$bp23_ag=$bp21 ."AG";
					# GC AND PLOY T
					my $count1=$bp20=~tr/G/G/;
					my $count2=$bp20=~tr/C/C/;
					my $gc = ($count1+$count2)/20;
					$escore=0;
					if(($gc<0.75)&&($gc>0.25)) {
						$escore= $escore+5;
					}
					if($bp20 =~ /TTTT/) {
						;
					}
					else {
						$escore= $escore+5;
					}
					
					printf OUTB ">$no_sgRNA\n$bp23\n";
					printf OUTD ">$no_sgRNA\n$bp23_ag\n";
					
					printf OUTA "$bp20\t$chromosome[$id_chr-2]\t$location\t+\t$escore\n"  ;
				}		
			}  ## while
			
		}
		if($id_chr>2) {
			$no_sgRNA_chr[$id_chr-2]= $no_sgRNA-$no_sgRNA_chr[$id_chr-3];
			printf OUTE "$chromosome[$id_chr-2]\t$no_sgRNA_chr[$id_chr-2]\n";	
		}	
		$seq="";  ##not sure
	
	}
	elsif ($_ =~ />/) {
		$seq = "";
	}
	else {
		my $temp=uc($_);		
		$seq=$seq.$temp;
	}
	
}

#### for the last chromatin
		if(length($seq) !=0) {
			my $no_gg=0;	
			
			my $flag=0;
			my $start_bp=21;
			# find sgRNA
			while($flag !=(-1)) {
				$flag = index($seq, "GG",$start_bp);
				my $start_23=0;
				if($flag !=(-1))  {$start_23=$flag-21;}
				$location=$start_23;
				$bp23 = substr ($seq,$start_23,23);
				$start_bp=$flag+1;
			
					$no_sgRNA ++;
					$bp20 = substr($bp23,0,20);
					$bp21 = substr($bp23,0,21); 
					$bp23_ag=$bp21 ."AG";
					
					my $count1=$bp20=~tr/G/G/;
					my $count2=$bp20=~tr/C/C/;
					my $gc = ($count1+$count2)/20;
					$escore=0;
					if(($gc<0.75)&&($gc>0.25)) {
						$escore= $escore+5;
					}
					if($bp20 =~ /TTTT/) {
						;
					}
					else {
						$escore= $escore+5;
					}
					
					printf OUTB ">$no_sgRNA\n$bp23\n";
					printf OUTD ">$no_sgRNA\n$bp23_ag\n";
					printf OUTA "$bp20\t$chromosome[$id_chr-1]\t$location\t-\t$escore\n";  ##location is not defined $chromosom
								
			}  ## while
			 $start_bp=0;
			$flag=0;
			$seq =$seq."MMMMMMMMMMMMMMMMMMMMMM";
			while($flag !=(-1)) {
				$flag = index($seq, "CC",$start_bp);
			
				my $start_23=0;
				if($flag !=(-1))  {$start_23=$flag;}
				$location=$start_23+3;
				$bp23 = substr ($seq,$start_23,23);
				
				$start_bp=$flag+1;
				if(substr($bp23,22,1) ne "M") {	
					$no_sgRNA ++;

					$bp20 = substr($bp23,3,20); 
					$bp21 = substr($bp23,2,21); 
					
					$bp20= reverse $bp20;
					$bp20=~tr/ACGT/TGCA/;
					$bp21= reverse $bp21;
					$bp21=~tr/ACGT/TGCA/;
					$bp23_ag=$bp21 ."AG";
					# GC AND PLOY T
					my $count1=$bp20=~tr/G/G/;
					my $count2=$bp20=~tr/C/C/;
					my $gc = ($count1+$count2)/20;
					$escore=0;
					if(($gc<0.75)&&($gc>0.25)) {
						$escore= $escore+5;
					}
					if($bp20 =~ /TTTT/) {
						;
					}
					else {
						$escore= $escore+5;
					}
					
					printf OUTB ">$no_sgRNA\n$bp23\n";
					printf OUTD ">$no_sgRNA\n$bp23_ag\n";
					
					printf OUTA "$bp20\t$chromosome[$id_chr-1]\t$location\t+\t$escore\n"  ;
				}		
			}  ## while
			
		}
		$no_sgRNA_chr[$id_chr-1]= $no_sgRNA-$no_sgRNA_chr[$id_chr-2];
		printf OUTE "$chromosome[$id_chr-1]\t$no_sgRNA_chr[$id_chr-1]\n";		

close INA;
close OUTA;
close OUTB;
close OUTC;
close OUTD;
close OUTE;


