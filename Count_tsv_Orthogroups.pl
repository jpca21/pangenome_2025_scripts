# !/usr/bin/perl
use strict;
use warnings;

#### Variables
my $file01 = $ARGV[0] || "new_Orthogroup_matrix.tsv" ;
my $outfile = $file01;
$outfile =~ s/\.\w+//i;

#### Abriendo el archivo (tipo Orthogroups.tsv)
my $line = 0;
my $scalar = 0;
open (A, $file01) or die;

#### Abriendo las salidas
open (OUT_CONT, ">$outfile\.Count.tsv") or die;
open (OUT_BIN,  ">$outfile\.Binnary.tsv") or die;

LOOOP:
while (my $l2 = <A>) {
    ++$line;
    $l2 =~ s/\n//g;
    $l2 =~ s/\r//g;

    #### La primera linea
    if ($line == 1) {
        print OUT_CONT "$l2\n";
        print OUT_BIN "$l2\n";
        my @data = split ("\t", $l2);
        $scalar = scalar(@data); ### Este es el maximo de columnas del archivo!
        next LOOOP;
    }

    #### Todas las demas
    my @data = split ("\t", $l2);
    my $curr_scalar = scalar(@data);
    #print scalar(@data); ### DIAG!
    #print "\n"; ### DIAG!
    my $diff = ($scalar - $curr_scalar) ; # Cuantas columnas le falta a la linea actual?
    #print $diff; ### DIAG!
    #print "\n"; ### DIAG!

    #### Rellenando
    my $curr_linea_index = 0;
    foreach my $dato (@data) {
       ++$curr_linea_index;
       if ($curr_linea_index == 1) { # la primera "celda es la descripcion"
           print OUT_CONT $dato;
           print OUT_BIN $dato;
       } elsif ($curr_linea_index != 1) { # Esto es para contar el contenido de la "celda" ->
           my $count = 0;
           my $binnr = 0;
           if (length($dato) < 2) { # si el largo del contenido de la "celda" es practicamente cero
               $count = 0;
               $binnr = 0;
               print OUT_CONT "\t$count";
               print OUT_BIN  "\t$binnr";
           } elsif (length($dato) > 2) { # Si la celda podria tener mas que cero elementos
               my @dss = split (",", $dato);
               my $cnt = scalar(@dss);
               $count = $cnt;
               $binnr = 1;
               print OUT_CONT "\t$count";
               print OUT_BIN  "\t$binnr";
           }
       }  # <-
    }

    #### Cuando ($diff > 0) hay que rellenar con zeros, por `$diff` veces
    for (my $i=0; $i < $diff; $i++) {
        print OUT_CONT "\t0";
        print OUT_BIN  "\t0";
    }

    #### Finalizando linea
    print OUT_BIN  "\n";
    print OUT_CONT "\n";

}

#### Final
close A;
close OUT_CONT;
close OUT_BIN;
exit;

