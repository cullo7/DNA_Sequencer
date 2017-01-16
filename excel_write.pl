#!/usr/bin/perl
use strict;
use warnings;

use Excel::Writer::XLSX;

my $workbook  = Excel::Writer::XLSX->new( 'output.xlsx' );
my $worksheet = $workbook->add_worksheet();

my $filename = 'output.txt';
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";

my $row = 0;
while (my $line = <$fh>) {
  chomp $line;
  my @values = split(',', $line);

  my $col = 0;
  foreach my $val (@values) {
	$worksheet->write($row,$col, $val );
    $col++;
  }
  $row++;
}

$workbook->close;
