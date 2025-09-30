# LAST VERSION: July 9, 2024

use strict;
use warnings;

#### other Modules
use Bio::SeqIO;

#### Inputs
my $OG_list = $ARGV[0] || "Orthogroups_SingleCopyOrthologues.txt";
my $extensive_list = $ARGV[1]  || "Orthogroups.txt";
my $extensive_fasta = $ARGV[2] || "ALL_FAA";
my $header_list = $ARGV[3] || "genome_to_id.tsv";

#### Configuring MAFFT
my $mafft_exe = "mafft";
my $options = "--maxiterate 1000 --localpair --op 5";

#### Obtaining OG list
print "### Obtaining OG list\n";
my @contentOG;
open (A, $OG_list) or die;
@contentOG = <A>;
close A;

#### Obtaining Extensive data (as a HASH)
print "### Obtaining Orthogroups to IDs data\n";
my %content_data;
open (B, $extensive_list) or die;
while (my $l = <B>) {
    $l =~ s/\n//i;
    $l =~ s/\r//i;
    # REGEXP
    $l =~ /(^\S+)\: (.+$)/ ;
    my $a = $1;
    my $b = $2;
    $content_data{$a} = $b;
}
close B;

#### Obtaining Headers list
print "### Obtaining Genome to IDs data\n";
my %contentHeaders;
open (C, $header_list) or die;
while (my $lx = <C>) {
    $lx =~ s/\n//i;
    $lx =~ s/\r//i;
    $lx =~ s/ //i;
    my ($colA, $colB) =  split ("\t", $lx);
    $contentHeaders{$colB} = $colA;
}
close C;

#### Creating HASH with fasta file content
print "### Obtaining FASTA data\n";
my %fasta = ();
my $inseq = Bio::SeqIO->new(-file   => $extensive_fasta,
                            -format => "Fasta", );
while (my $seq = $inseq->next_seq) {
  my $aa = $seq->id;
  my $bb = $seq->seq;
  (chop $bb) if ($bb =~ /\*$/);
  ### in case of format "ACC|LocusTag"
  #my @x = split ("\|", $aa);
  #my $locus =$x[1];
  #$fasta{$locus} = $bb;
  $fasta{$aa} = $bb; ### OLD

}

#### TEST FASTA DATA
# foreach my $qwq ( keys(%fasta)) {
#   my $qeq = $fasta{$qwq};
#   print "$qwq\t$qeq\n";
# }
# print "\n";

#### TEST Orthogroup.txt data
# foreach my $cwq ( keys(%content_data)) {
#   my $ceq = $content_data{$cwq};
#   print "$cwq\t$ceq\n";
# }
# print "\n";


#### TEST genometoid.txt data
# foreach my $cww ( keys(%contentHeaders)) {
#   my $cew = $contentHeaders{$cww};
#   print "$cww\t$cew\n";
# }
# print "\n";

#########################################################################################
#### VARIANT
#### NOW, we will start the loop to create the Concatenator
#### 1- for each select OG
####  2- for each HEADER
####    3- search seqIds in the extensive list
####      4- use found seqIds in the fasta hash
####        5- align current OG profile
####          6- Open ALN, store and concatentate seq string in ordered form

#### hash for writing new concatenated Fasta
my %new_fasta_conc = ();
my $conc_seq = "";

A01:
foreach my $curr_OG (@contentOG) {
    $curr_OG =~ s/\n//;
    $curr_OG =~ s/\r//;
    $curr_OG =~ s/ //g;

    print "Aligning Block: $curr_OG...\n";

    open (X, ">$curr_OG.txt") or die;

    #### Searching inside the selected OGs
    if (exists $content_data{$curr_OG}) { #### OG!
        my $d = $content_data{$curr_OG};
        my @data = split (" ", $d); ### All seqs in the OG

        A02: ### Focusing in each seq
        foreach my $e (@data) { ### For each id
            #### Searching header (index substr in str)
            my $ff = $contentHeaders{$e};
            $ff =~ s/ //g;
            if (exists $fasta{$e}) {
                my $sqq = $fasta{$e};
                print X "\>$ff\n$sqq\n"; ### DIAG
            }
        }
    }
    close X;

  #### Creating ALIGNMENT TEMPORAL (MAFFT must be executable as environment)
  system ("$mafft_exe --quiet $options $curr_OG.txt > $curr_OG.aln");

  #### Storing aln and include it into the current HASH
  my $inseq2 = Bio::SeqIO->new(-file   => "$curr_OG.aln",
                              -format => "Fasta", );
  while (my $seq2 = $inseq2->next_seq) {
    my $aaa = $seq2->id;
    my $bbb = $seq2->seq;
# my
    if(exists $new_fasta_conc{$aaa}) { ## if already present
      my $ccc = $new_fasta_conc{$aaa};
      $new_fasta_conc{$aaa} = $ccc . $bbb;
    } else { ## if it is the first
      $new_fasta_conc{$aaa} = $bbb;
    } #
  }
  #die; #DIAG
  #### Delete temporal files
  unlink ("$curr_OG.txt");
  unlink ("$curr_OG.aln");

  #die; ### DIAG


} ## END A01


#########################################################################################
### Final parsing
my $OUTPUT = $ARGV[4] || "OUTPUT_CONC.ALN.fasta";
open (OUT, ">$OUTPUT") or die ;

print "Assembling...\n";
while ((my $key3, my $value3) = each (%new_fasta_conc)) {
    print  OUT "\>$key3\n" . $new_fasta_conc{$key3} .  "\n";
}

close OUT;
print "Done! output available in $OUTPUT\n";
#### END
exit;


#########################################################################################
### Annex: functions
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
#########################################################################################
