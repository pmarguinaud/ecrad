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

sub addArgLACC
{
  my $pu = shift;

  my ($acc_routine) = &F ('./routine-openacc', $pu);

  return if ($acc_routine);

  my $stmt = $pu->firstChild;

  my ($arglist) = &F ('./dummy-arg-LT', $stmt);

  if ($arglist)
    {
      my @arg = &F ('./arg-N', $arglist, 1);
      return if (grep { lc ($_) eq 'lacc' } @arg);
      $arglist->appendChild (&t (', '));
    }
  else
    {
      $arglist = &n ('<dummy-arg-LT/>');
      $stmt->appendChild ($_) for (&t ('('), $arglist, &t (')'));
    }

  $arglist->appendChild (&n ('<arg-N><N><n>lacc</n></N></arg-N>'));

  my @decl = &F ('./T-decl-stmt', $pu);

  my $last = $decl[-1];

  unless ($last)
     {
       $last = $stmt;
     }

  my $ind = "\n" . (' ' x &Fxtran::getIndent ($last));


  for (&s ("logical, intent (in) :: lacc"), &t ($ind))
    {
      $last->parentNode->insertAfter ($_, $last);
    }
}

    
sub process
{
  my $d = shift;

  my @pu = &F ('.//program-unit', $d);
  
  for my $pu (@pu)
    {
      my $stmt = $pu->firstChild;
      if ($stmt->nodeName =~ m/^(?:subroutine|function)-stmt$/o)
        {
          &addArgLACC ($pu);
        }
    }
  
  
  my @skip_call = qw (fast_adding_ica_lw calc_fluxes_no_scattering_lw adding_ica_sw cloud_generator_acc indexed_sum
                      calc_ice_optics_fu_sw calc_ice_optics_fu_lw calc_liq_optics_slingo calc_liq_optics_lindner_li
                      calc_liq_optics_socrates calc_lw_derivatives_ica modify_lw_derivatives_ica calc_two_stream_gammas_lw
                      calc_two_stream_gammas_sw calc_reflectance_transmittance_lw calc_ref_trans_lw
                      calc_no_scattering_transmittance_lw calc_reflectance_transmittance_sw calc_ref_trans_sw);

  my %skip_call = map { (lc ($_), 1) } @skip_call;

  my @call = &F ('.//call-stmt', $d);

  my $getlaccval = sub
  {
    my $node = shift;
    my @pu = &F ('ancestor::program-unit', $node);
    my $pu = pop (@pu);
    
    my $laccval = 'lacc';

    my $stmt = $pu->firstChild;
    if ($stmt->nodeName eq 'subroutine-stmt')
      {
        my ($N) = &F ('./subroutine-N', $stmt, 1);
        $laccval = '.false.' if (lc ($N) eq 'read_input_gpu');
      }

    return &e ($laccval);
  };
  
  for my $call (@call)
    {
      my $laccval = $getlaccval->($call);

      my ($proc) = &F ('./procedure-designator', $call, 1);
      
      $proc =~ s/.*%//o;

      next unless ($proc=~ s/_gpu$//io);

      next if ($skip_call{lc ($proc)});
  
      my ($argspec) = &F ('./arg-spec', $call);
  
      if ($argspec)
        {
          my @arg = &F ('./arg', $argspec, 1);
          next if (grep { (lc ($_) eq 'lacc') || (lc ($_) eq 'lacc=lacc')  || (lc ($_) eq 'lacc=llacc') } @arg);

          if (grep { lc ($_) eq 'lacc=.true.' } @arg)
            {
              my ($arg) = grep { lc ($_->textContent) =~ m/^lacc=/o } &F ('./arg', $argspec);
              $arg->replaceNode (&n ("<arg><arg-N><k>lacc</k></arg-N>=$laccval</arg>"));
              next;
            }
          if (grep { lc ($_) eq 'lacc=.false.' } @arg)
            {
              next;
            }

          $argspec->appendChild (&t (', ')) if (@arg);
        }
      else
        {
          $argspec = &n ("<arg-spec/>");
          $call->appendChild ($_) for (&t ('('), $argspec, &t (')'));
        }
  
      $argspec->appendChild (&n ("<arg><arg-N><k>lacc</k></arg-N>=$laccval</arg>"));
    }
  
  my @expr = &F ('.//named-E[./R-LT/parens-R]', $d);
  
  my @skip_func = qw (calc_rh_index beta2alpha initialize_acc uniform_distribution_acc);
  my %skip_func = map { (lc ($_), 1) } @skip_func;

  for my $expr (@expr)
    {
      my $laccval = $getlaccval->($expr);

      my ($N) = &F ('./N', $expr, 1);
      next unless ($N =~ s/_gpu$//io);
      next if ($skip_func{lc ($N)});
      my ($elt) = &F ('./R-LT/parens-R/element-LT', $expr);
      $elt->appendChild (&t (', '));
      $elt->appendChild (&n ("$laccval"));
    }

  my @meth = qw (find min_wavenumber max_wavenumber out_of_physical_bounds);

  for my $meth (@meth)
    {
      my @n = (&F ('.//component-R/ct/text()[translate(string(.),"abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ")="?"]', uc ($meth . '_gpu'), $d));
      for my $n (@n)
        {
          my $laccval = $getlaccval->($n);

          my $c = $n->parentNode->parentNode->parentNode;
          my @r = &F ('./ANY-R', $c);
          my $p = $r[-1];
          my ($elt) = &F ('./element-LT', $p);
          $elt->appendChild (&t (', ')) if ($elt->childNodes ());
          $elt->appendChild (&n ("$laccval"));
        }
    }

}


my $f = shift;

if ($f =~ m/\.F90$/o)
  {
    my $F90 = $f;
    
    $F90 =~ s/\.F90$//o;
    $F90 =~ s/_gpu//o;
    $F90 = "$F90\_gpu.F90";

    
    my $d = &Fxtran::parse (location => $F90, fopts => [qw (-construct-tag -no-include -line-length 500 -openacc)]);

    &process ($d);

#   print $d->textContent;
    'FileHandle'->new (">$F90")->print ($d->textContent);
  }
elsif ($f =~ m/\.intfb\.h$/o)
  {
    my $intfb = $f;

    $intfb =~ s/\.intfb\.h$//o;
    $intfb =~ s/_gpu//o;
    $intfb = "$intfb\_gpu.intfb.h";

    my ($d) = &fxtran::parse (fragment => &slurp ($intfb));
    &process ($d);

    'FileHandle'->new (">$intfb")->print ($d->textContent);
  }


