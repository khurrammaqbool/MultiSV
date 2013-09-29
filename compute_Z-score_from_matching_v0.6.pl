#!/usr/bin/perl -w
use strict;

my $F3_readfile = $ARGV[0];
my $R3_readfile = $ARGV[1];
my $Z_score_file = $ARGV[2];
open F3_FILE, $F3_readfile or die "$F3_readfile $!";
open R3_FILE, $R3_readfile or die "$R3_readfile $!";
open (Z_SCORE_FILE,">$Z_score_file") or die "$Z_score_file $!";

my $chromosome_in_F3_file;
my $head_of_F3_file = "head_of_F3_file";
my $cmd = "head $F3_readfile > $head_of_F3_file";
system("head $F3_readfile > $head_of_F3_file");
open F3_FILE, $head_of_F3_file or die "$head_of_F3_file $!";
while(<F3_FILE>)
{
    chomp;s/,/\t/g;s/:\(/\t\(/g;s/>//;
    my @header = split;my @line = ();my $read_start;my $read_end;
    if(defined $header[1])
      {$header[1] =~ s/\./\t/;$header[1] =~ s/_/\t/;
           if ($header[1] =~ /-/){$header[1] =~ s/\-/\-\t/;}
           else  {$header[1] =~  s/\t/\t+\t/ ;}
     @line = split '\t', $header[1];
       }
      my $orientation = $line[1];my $position = $line[2];
      if(defined $position)
        {$chromosome_in_F3_file = $line[0];}
}
print $chromosome_in_F3_file,"\n";










my %average_genome_window_coverage;
my %chromosome_size;
$chromosome_size{1} = 333917662;$chromosome_size{2} = 169113706 ;$chromosome_size{3} = 157499821;$chromosome_size{4} = 156303021;
$chromosome_size{5} = 122295561 ;$chromosome_size{6} = 160888347;$chromosome_size{7} = 399598981 ;$chromosome_size{8} = 157303522;
$chromosome_size{9} = 164618141;$chromosome_size{10} = 86567162;$chromosome_size{11} = 92581275 ;$chromosome_size{12} = 69842989 ; 
$chromosome_size{13} = 161749637;$chromosome_size{14} = 175947320;$chromosome_size{15} = 164292966;$chromosome_size{16} = 92518118 ;
$chromosome_size{17} = 76063987;$chromosome_size{18} = 64311137;$chromosome_size{19} = 143508085;$chromosome_size{20} =  1061286;
$chromosome_size{MT} =  16855;

my %position_coverage=();
my $position_coverage_string_1 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
my $position_coverage_string_2 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
my $position_coverage_string_3 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
my $position_coverage_string_4 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
my $position_coverage_string_5 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
my $position_coverage_string_6 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
#my $position_coverage_string_7 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
#my $position_coverage_string_8 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
#my $position_coverage_string_9 = chr(0)x$chromosome_size{$chromosome_in_F3_file};
#my $position_coverage_string_10 = chr(0)x$chromosome_size{$chromosome_in_F3_file};


my $window_size = 1000;
my $genome_window = 0;
my %genome_window_coverage=();
my $number_of_reads_covered;


while(<F3_FILE>)
{
    chomp; 
    s/,/\t/g;
    s/:\(/\t\(/g;
    s/>//;
    my @temp = split;
    my @line = ();
    my $read_start;
    my $read_end;
    
    if(defined $temp[1])
      {
        $temp[1] =~ s/\./\t/;
        $temp[1] =~ s/_/\t/;
           if ($temp[1] =~ /-/)
              {$temp[1] =~ s/\-/\-\t/;}
           else    
              {$temp[1] =~  s/\t/\t+\t/ ;}
          	@line = split '\t', $temp[1];
          	#print join("\t", $line[2], "\n") ;
                 #print join("\t", @line, "\n") if scalar @temp == 3;
       }
                
      my $orientation = $line[1];
      my $position = $line[2];
            
      if(defined $position)
        {
           if ($orientation eq '+')
             {
				($read_start,$read_end) = ($position+1,$position+50 -1);
              }
          elsif ($orientation eq '-') 
          	 {
    			  $read_start = $position - 50 + 3;
   				  $read_end = $position + 1;
		  	 }                       
            
		for (my $i=$read_start; $i < $read_end; $i++)
		{   my $position_cov;
				if (ord(substr($position_coverage_string_1,($i-1),1)) < 255)
				{ 
					$position_cov = ord(substr($position_coverage_string_1, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_1, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_2,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_2, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_2, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_3,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_3, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_3, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_4,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_4, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_4, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_5,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_5, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_5, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_6,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_6, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_6, ($i-1),1)) = chr($position_cov);
				}

				
			}
		}           
}
print "File F3 Complete\n";
while(<R3_FILE>)
{
    chomp; 
    s/,/\t/g;
    s/:\(/\t\(/g;
    s/>//;
    my @temp = split;
    my @line = ();
    my $read_start;
    my $read_end;
    if(defined $temp[1])
      {
        $temp[1] =~ s/\./\t/;
        $temp[1] =~ s/_/\t/;
           if ($temp[1] =~ /-/)
              {$temp[1] =~ s/\-/\-\t/;}
           else    
              {$temp[1] =~  s/\t/\t+\t/ ;}
          	@line = split '\t', $temp[1];
       }
                
      my $orientation = $line[1];
      my $position = $line[2];
      if(defined $position)
        {
          if ($orientation eq '+')
             {
	    	($read_start,$read_end) = ($position+1,$position+50 -1);
              }
          elsif ($orientation eq '-') 
          	 {
    		   $read_start = $position - 50 + 3;
   	       	  $read_end = $position + 1;
		  }                       
            
		for (my $i=$read_start; $i < $read_end; $i++)
		{ my  $position_cov;
				if (ord(substr($position_coverage_string_1,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_1, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_1, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_2,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_2, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_2, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_3,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_3, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_3, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_4,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_4, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_4, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_5,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_5, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_5, ($i-1),1)) = chr($position_cov);
				}
				elsif (ord(substr($position_coverage_string_6,($i-1),1)) < 255)
				{
					$position_cov = ord(substr($position_coverage_string_6, ($i-1),1)) + ord(chr(1));
					(substr($position_coverage_string_6, ($i-1),1)) = chr($position_cov);
				}
      	       	}
		
     }
}   

print "File R3 Complete\n";	
close F3_FILE; close R3_FILE;

	for (my $i=1; $i < $chromosome_size{$chromosome_in_F3_file}; $i=$i + $window_size)
#for ($i=1; $i < 1000000; $i=$i + $window_size)
{  
	for (my $position = $i; $position < $i + $window_size and $position < $chromosome_size{$chromosome_in_F3_file} ; $position++)
		{
			$genome_window_coverage{$genome_window} += ord(substr($position_coverage_string_1, ($position-1),1)) + 
			 ord(substr($position_coverage_string_2, ($position-1),1)) + ord(substr($position_coverage_string_3, ($position-1),1))
			 +  ord(substr($position_coverage_string_4, ($position-1),1)) 
			 +  ord(substr($position_coverage_string_5, ($position-1),1))
			 + ord(substr($position_coverage_string_6, ($position-1),1))
			 ;
			 $average_genome_window_coverage{$genome_window} = ($genome_window_coverage{$genome_window} / $window_size);
				
		}
			$genome_window++;			
}


print "Calculation of average coverage for each window  Complete\n";

my $combined_sum;
for (my $i = 0; $i < $genome_window; $i++)
{
	 $combined_sum += $average_genome_window_coverage{$i};
	
}
	my $combined_average = $combined_sum / $genome_window;
my $squared_deviation_from_mean;
for (my $i = 0; $i < $genome_window; $i++)
{
	 $squared_deviation_from_mean +=  ( $combined_average - $average_genome_window_coverage{$i} )**2;
}
my $variance = ($squared_deviation_from_mean / ($genome_window -1 ));

my $standard_deviation = sqrt ($variance);
my $z_score;
for (my $i = 0; $i < $genome_window; $i++)
{
	if ($average_genome_window_coverage{$i} != 0)
       
	 
	{$z_score = sprintf "%.2f",(($average_genome_window_coverage{$i} - $combined_average) / $standard_deviation);}
	else 
	{$z_score = "NA";}
	print Z_SCORE_FILE $chromosome_in_F3_file,"\t",$i,"\t",$z_score,"\n";
}

print "Z_score calculation  Complete\n";


=pod


my $R_script;
my $R_script_continu;
my $scriptfile = "RPlotScript.R";

$R_script = "genome_window <- read.table('$outfile', sep = \"\\t\")\n";
$R_script.="library(plyr)\n";
$R_script.="splitdel <- splitter_d(genome_window, .(genome_window\$V1))\n";
$R_script.="attach(mtcars)\n";
$R_script.="par(mfrow=c(5,6))\n";

open OUT, ">$scriptfile" or die "$scriptfile $!";
print OUT $R_script;


for ($i=1,$p=1; $i < $genome_window +1 ; $i++)
{ 
	if ($p < $genome_window and ($p/$i==1 or $i==1))
	{  
	
	$p=$p+30;
#	$R_script_continu=<<END;
#	dev.off()
#	END
	$R_script_continu="dev.off()\n";
	print OUT $R_script_continu;
	
	$R_script_continu="png(file=\"png/$outfile$p.png\")\n";
	print OUT $R_script_continu;
	$R_script_continu="par(mfrow=c(5,6))\n";
#	$R_script_continu="par(mar=c(2,2,0,0))\n";
	print OUT $R_script_continu;
	}
	
	
#	$R_script_continu="barplot(splitdel[[$i]]\$V2,main='del$1',ylab=\"Cov\",xlab=\"bp_position$i\")\n";
	$R_script_continu="barplot(splitdel[[$i]]\$V2,ylab=\"Cov\",xlab=\"$i\")\n";
	print OUT $R_script_continu;
}	
close OUT;



my $cmd = 'R CMD BATCH RPlotScript.R';
system($cmd);



=cut
