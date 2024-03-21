#!/usr/bin/perl -w

use strict;

use FileHandle;
use Data::Dumper;
use Getopt::Long;
use File::stat;
use File::Path;
use File::Copy;
use File::Basename;
use FindBin qw ($Bin);
use lib "$Bin/../lib";

use lib "$ENV{HOME}/fxtran-acdc/lib";

use Common;
use Fxtran;

sub slurp
{
  my $f = shift;
  my $text = do { my $fh = 'FileHandle'->new ("<$f"); local $/ = undef; <$fh> };
  return $text;
}

for my $F90 (@ARGV)
  {
    my $d = &Fxtran::parse (location => $F90, fopts => [qw (-construct-tag -no-include -line-length 500)]);

    my @func = &F ('.//function-stmt', $d);

    for my $func (@func)
      {
        my ($N) = &F ('./ANY-N', $func);
        print $N->textContent, "\n";
      }
   }


