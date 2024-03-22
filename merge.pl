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

my @suff = qw (_gpu _cpu);
my $DIR = '/home/sor/ecrad/cpu';

sub slurp
{
  my $f = shift;
  my $text = do { my $fh = 'FileHandle'->new ("<$f"); local $/ = undef; <$fh> };
  return $text;
}

sub insertMethods
{
  my ($t, $suff, @p) = @_;

  my $end = $t->lastChild;

  my $indent = &Fxtran::getIndent ($end);

  for my $p (@p)
    {
      my @N = (&F ('./procedure-N-LT/rename/use-N/n/text()', $p),
               &F ('./procedure-N-LT/rename/N/n/text()', $p));
      for my $N (@N)
        {
          $N->setData ($N->data . uc ($suff));
        }

      $t->insertBefore ($p, $end);
      $t->insertBefore (&t ("\n\n" . (' ' x $indent)), $end);
    }
}

sub insertGenerics
{
  my ($t, $suff, @p) = @_;

  my $end = $t->lastChild;

  my $indent = &Fxtran::getIndent ($end);

  for my $p (@p)
    {
      my @N = (&F ('.//generic-spec/N/n/text()', $p),
               &F ('.//binding-N/N/n/text()', $p));
      for my $N (@N)
        {
          $N->setData ($N->data . uc ($suff));
        }

      $t->insertBefore ($p, $end);
      $t->insertBefore (&t ("\n\n" . (' ' x $indent)), $end);
    }
}

sub mergeType
{
  my ($t0, $t1, $t2) = @_;

  my @c0 = &F ('./component-decl-stmt', $t0);
  my @c1 = &F ('./component-decl-stmt', $t1);
  my @c2 = &F ('./component-decl-stmt', $t2);

  die unless (scalar (@c0) == scalar (@c1));
  die unless (scalar (@c0) == scalar (@c2));
   
  for my $i (0 .. $#c0)
    {
      die unless ($c0[$i]->textContent eq $c1[$i]->textContent);
      die unless ($c0[$i]->textContent eq $c2[$i]->textContent);
    }

  my @p0 = &F ('./procedure-stmt', $t0);
  my @p1 = &F ('./procedure-stmt', $t1);
  my @p2 = &F ('./procedure-stmt', $t2);

  $_->unbindNode ()
    for (@p0);

  &insertMethods ($t0, $suff[0], @p1);
  &insertMethods ($t0, $suff[1], @p2);

  my @g0 = &F ('./generic-stmt', $t0);
  my @g1 = &F ('./generic-stmt', $t1);
  my @g2 = &F ('./generic-stmt', $t2);

  $_->unbindNode ()
    for (@g0);

  &insertGenerics ($t0, $suff[0], @g1);
  &insertGenerics ($t0, $suff[1], @g2);
}

sub mergeF90
{
  my $F90 = shift;
  
  my @F90 = (map { my $f90 = $F90; my $suff = $_; $f90 =~ s/\.F90$/$suff.F90/o; $f90 } @suff);
  
  unshift (@F90, $F90[0]);
  
  my @d = map 
  { 
    my $F90 = $_; 
    &Fxtran::parse (location => $F90, fopts => [qw (-construct-tag -no-include -line-length 500)]) 
  } @F90;
  
  my ($pu0) = &F ('.//program-unit[1]', $d[0]);
  my ($pu1) = &F ('.//program-unit[1]', $d[1]);
  my ($pu2) = &F ('.//program-unit[1]', $d[2]);
  
  
  my $first = $pu0->firstChild;
  my $last = $pu0->lastChild;
  
  if ($first->nodeName eq 'module-stmt')
    {
      my @pu0 = &F ('./program-unit', $pu0);
      my @pu1 = &F ('./program-unit', $pu1);
      my @pu2 = &F ('./program-unit', $pu2);
      
      for (@pu0)
        {
          $_->unbindNode ();
        }
      
      my $end = $pu0->lastChild;
      
      
      for my $pu (@pu1, @pu2)
        {
          my $indent = &Fxtran::getIndent ($pu);
          $pu0->insertBefore (&t (' ' x $indent), $end);
          $pu0->insertBefore ($pu, $end);
          $pu0->insertBefore (&t ("\n\n"), $end);
        }
      
      
      my @t0 = &F ('./T-construct', $pu0);
      my @t1 = &F ('./T-construct', $pu1);
      my @t2 = &F ('./T-construct', $pu2);
      
      die unless (scalar (@t0) == scalar (@t1));
      die unless (scalar (@t0) == scalar (@t2));
      
      for my $i (0 .. $#t0)
        {
          &mergeType ($t0[$i], $t1[$i], $t2[$i]);
        }

      for my $stmt ($first, $last)
        {
          my ($N) = &F ('./ANY-N/N/n/text()', $stmt);
          if ($N)
            {
              (my $tt = $N->data) =~ s/\_(CPU|GPU)$//o;
              $N->setData ($tt);
            }
        }

      
      &mkpath (&dirname ("$DIR/$F90"));
      'FileHandle'->new (">$DIR/$F90")->print ($d[0]->textContent);
    }
  elsif ($first->nodeName eq 'program-stmt')
    {
      &mkpath (&dirname ("$DIR/$F90"));
      
      my $n;

      for my $pu ($pu1, $pu2)
        {
          my ($first, $last) = ($pu->firstChild, $pu->lastChild);

          my ($N) = &F ('./ANY-N', $first, 1);

          $first->replaceNode (&n ("<subroutine-stmt>SUBROUTINE <subroutine-N><N><n>$N</n></N></subroutine-N></subroutine-stmt>"));
          $last->replaceNode (&n ("<end-subroutine-stmt>END SUBROUTINE</end-subroutine-stmt>"));

          ($n = $N) =~ s/_(?:CPU|GPU)$//o;
        }
     
      my $fh = 'FileHandle'->new (">$DIR/$F90");

      $fh->print (<< "EOF");
PROGRAM $n

IMPLICIT NONE

CHARACTER*16 :: CLARCH

CALL GETENV ('ARCH', CLARCH)

SELECT CASE (CLARCH)
  CASE ('CPU') 
    CALL ${n}_CPU
  CASE ('GPU') 
    CALL ${n}_GPU
  CASE DEFAULT
    ERROR STOP "UNEXPECTED ARCH"
END SELECT

CONTAINS
EOF

      $fh->print ("\n\n");
      $fh->print ($pu1->textContent);
      $fh->print ("\n\n");
      $fh->print ($pu2->textContent);
      $fh->print ("\n\n");

      $fh->print ("END PROGRAM $n\n");

      $fh->close ();
    }
  elsif ($first->nodeName eq 'subroutine-stmt')
    {
      &mkpath (&dirname ("$DIR/$F90"));
      'FileHandle'->new (">$DIR/$F90")->print (join ("\n", $d[1]->textContent, "", $d[2]->textContent));
    }
  else
    {
      die;
    }
}

sub mergeIntfb
{
  my $intfb = shift;
  
  my @intfb = (map { my $f = $intfb; my $suff = $_; $f =~ s/\.intfb.h$/$suff.intfb.h/o; $f } @suff);
  
  'FileHandle'->new (">$DIR/$intfb")->print (join ("\n", map { &slurp ($_) } @intfb));
}

my ($file) = @ARGV;


if ($file =~ m/\.F90$/o)
  {
    &mergeF90 ($file);
  }
elsif ($file =~ m/\.intfb.h$/o)
  {
    &mergeIntfb ($file);
  }
else
  {
    die;
  }
