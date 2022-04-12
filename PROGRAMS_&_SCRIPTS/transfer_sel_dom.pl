#!/usr/local/bin/perl
use Getopt::Long;
use File::Find;
use strict;

my (%options,@TABLE, @vector, $x);

#options
GetOptions(
     'i=s'  => \$options{'i'},
     't=s'  => \$options{'t'},
     'o=s'  => \$options{'o'},
     'p=s'  => \$options{'p'}
);

#Usage
if($options{'i'} eq undef || $options{'t'} eq undef || $options{'o'} eq undef  || $options{'p'} eq undef) {
   print "\nUsage: \nperl transfer_sel_dom.pl -i (slim outfile)  -t (table sel and dom) -o (output file) -p (1/0 to include phenotype effect)\n\n";
   exit;
}

#read the table file
open(INPUTT,"$options{'t'}") or die;
$_=<INPUTT>; #skip first row
while(<INPUTT>) {
   $_ =~ /(\d+)\s(-?\d+\.?\d+)\s(-?\d+\.?\d+)\s(-?\d+\.?\d+)/; 
   @vector = ($1,$2,$3,$4);
   push (@TABLE,@vector); #keep all in the array TABLE
}
close(INPUTT);

#Read slim file and print new slim file with selection and dominance from TABLE
open(INPUTS,"$options{'i'}") or die;
open(OUTPUT,">$options{'o'}") or die;
while(<INPUTS>) {
   if($_=~ /(\d+)\s(\d+)\sm2\s(\d+)\s0\s0\sp1\s(\d+)\s(\d+)/ ) {
      for($x=0; $x<$#TABLE/4; $x++) {
         if($3==$TABLE[$x*4]) {
            if($options{'p'} eq '0') {
               print OUTPUT $1." ".$2." "."m2"." ".$3." ".$TABLE[$x*4+1]." ".$TABLE[$x*4+2]." "."p1"." ".$4." ".$5."\n";
            }else {
               print OUTPUT $1." ".$2." "."m2"." ".$3." ".$TABLE[$x*4+1]." ".$TABLE[$x*4+2]." "."p1"." ".$4." ".$5." ".$TABLE[$x*4+3]."\n";
            }
            $x=$#TABLE; #end and break
         }
      }
   }
   else {
      print OUTPUT $_;
   }
}
close(INPUTS);
close(OUTPUT);



