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

sub addSuffix
{
  my ($d, $suff, $openacc) = @_;

  my @pu = &F ('.//program-unit', $d);
  
  for my $pu (@pu)
    {
      if ($pu->firstChild->nodeName eq 'function-stmt')
        {
          my ($N) = &F ('./ANY-N/N/n/text()', $pu->firstChild, 2);
          if (my ($EN) = &F ('.//EN-N/N/n/text()[string(.)="?"]', $N, $pu))
            {
              $EN->setData ($EN->data . uc ($suff));
            }
          my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $N, $pu);
          for my $expr (@expr)
            {
              $expr->setData ($expr->data . uc ($suff));
            }
        }
      for my $f ($pu->firstChild, $pu->lastChild)
        {
          my ($N) = &F ('./ANY-N/N/n/text()', $f);
          $N->setData ($N->textContent . uc ($suff))
            if ($N);
        }
    }
  
  for (&F ('.//unseen', $d))
    {
      $_->unbindNode ();
    }
  
  unless ($openacc)
    {
      for (&F ('.//C[starts-with (string (.),"!$ACC")]', $d))
        {
          $_->unbindNode ();
        }
    }
  
  my @call = &F ('.//call-stmt', $d);
  
  for my $call (@call)
    {
      my ($proc) = &F ('./procedure-designator/named-E', $call);
      if (my @r = &F ('./R-LT/ANY-R', $proc))
        {
          my ($N) = &F ('./N', $proc, 1);
          
          unless ($N =~ m/^(?:FILE|OUT_FILE|driver_config)$/io)
            {
               my @skip = qw (load_netcdf);
               my %skip = map { (uc ($_), 1) } @skip;
               
               my $r = $r[-1];
               my ($t) = &F ('./ct/text()', $r);

               unless ($skip{uc ($t->data)})
                 {
                   $t->setData ($t->data . uc ($suff));
                 }
            }
        }
      else
        {
          my @skip = qw (DR_HOOK ABOR1 nvtxEndRange nvtxStartRange random_seed GETENV acc_attach get_command_argument radiation_abort 
                         uniform_distribution delta_eddington delta_eddington_extensive delta_eddington_extensive_vec
                         delta_eddington_scat_od revert_delta_eddington surrtab initialize_random_numbers surrtftr modify_wv_continuum
                         diag_mat_right_divide_3 expm expm_lw expm_sw fast_expm_exchange_2 fast_expm_exchange_3 fast_expm_exchange_3_orig
                         identity_minus_mat_x_mat identity_minus_mat_x_mat_3_lw identity_minus_mat_x_mat_3_sw lu_factorization lu_factorization_lw
                         lu_substitution lu_substitution_lw_inplace mat_square_sw mat_x_mat mat_x_mat_3_lw mat_x_mat_3_sw mat_x_mat_6 mat_x_mat_sw
                         mat_x_mat_sw_repeats mat_x_singlemat mat_x_singlemat_3_sw mat_x_vec mat_x_vec_3 mat_x_vec_3_lw
                         mat_x_vec_3_sw repeated_square repeated_square_sw_9 singlemat_x_mat singlemat_x_mat_3_sw singlemat_x_vec singlemat_x_vec_lw
                         singlemat_x_vec_sw solve_mat solve_mat_2 solve_mat_3 solve_mat_3_inplace solve_mat_3_lw solve_mat_3_sw solve_mat_inplace
                         solve_mat_n solve_mat_sw solve_vec solve_vec_2 solve_vec_3 solve_vec_3_lw solve_vec_3_lw_inplace solve_vec_3_ng
                         solve_vec_3_sw sparse_x_dense fast_expm_exchange random_number load_and_init INIT_SPECTRAL_PLANCK
                         srtm_kgb16 srtm_kgb17 srtm_kgb18 srtm_kgb19 srtm_kgb20 srtm_kgb21 srtm_kgb22 
                         srtm_kgb23 srtm_kgb24 srtm_kgb25 srtm_kgb26 srtm_kgb27 srtm_kgb28 srtm_kgb29);



          my %skip = map { (uc ($_), 1) } @skip;
  
          my ($t) = &F ('./N/n/text()', $proc);

          
          my $skip = $skip{uc ($t->data)};
 
          $skip ||= $t->data =~ m/^RRTM_KGB\d+/io;
  
          unless ($skip)
            {
              $t->setData ($t->data . uc ($suff));
            }
  
        }
    }

  # <C>!$omp declare simd(sample_from_pdf_simd) unifor

  my @OMP = &F ('.//C[starts-with(string(.),"!$omp")]/text()', $d);

  for my $OMP (@OMP)
    {
      $OMP->unbindNode
        if ($OMP->textContent =~ m/simd|linear/io);
    }

  # named-E><N><n>uniform_distribution_acc</n></N> <R-LT><parens-R
  
  my @func = qw (uniform_distribution_acc calc_beta_overlap_matrix calc_alpha_overlap_matrix
  calc_beta_overlap_matrix_dp calc_alpha_overlap_matrix_dp calc_planck_function_wavenumber
  get_line calc_rh_index string_loc_in_array lower_case DRY_AEROSOL_MASS_EXTINCTION planck_function indrad
  beta2alpha cloud_cover);


  
  for my $func (@func)
    {
      my @n = &F ('.//named-E[string(N)="?"]/N/n/text()', $func, $d);
      for my $n (@n)
        {
          $n->setData ($n->data . uc ($suff));
        }
    }

  my @meth = qw (find min_wavenumber max_wavenumber calc_rh_index out_of_physical_bounds);

  for my $meth (@meth)
    {
      my @n = (&F ('.//component-R/ct/text()[string(.)="?"]', lc ($meth), $d),
               &F ('.//component-R/ct/text()[string(.)="?"]', uc ($meth), $d));
      for my $n (@n)
        {
          $n->setData ($n->data . uc ($suff));
        }
    }

}

  
sub splitF90
{
  my $F90 = shift;

  for my $openacc (0, 1)
    {
      my $suff = $openacc ? '_gpu' : '_cpu';
  
      my $d = &Fxtran::parse (location => $F90, fopts => [qw (-construct-tag -no-include -line-length 500), 
                              $openacc ? qw (-D_OPENACC): ()]);

      $_->unbindNode ()
        for (&F ('.//cpp-section', $d));

      for my $cpp (&F ('.//cpp[@ref]', $d))
        {
          my $p = $cpp->parentNode;
          for ($cpp->childNodes ())
            {
              $p->insertBefore ($_, $cpp);
            }
          $cpp->unbindNode ();
        }

      my @use = &F ('.//use-stmt', $d);

      for my $use (@use)
        {
          my ($N) = &F ('./module-N', $use, 1);
          next if ($N =~ m/^yoe/io);
          my @c = $use->childNodes ();
          shift (@c) for (1 .. 2);
          $_->unbindNode for (@c);
        }

      my @public = &F ('.//public-stmt', $d);
      my @private = &F ('.//private-stmt', $d);

      $_->unbindNode 
        for (@public, @private);
  
      &addSuffix ($d, $suff, $openacc);

      (my $f90 = $F90) =~ s/\.F90$/$suff.F90/;
      
      'FileHandle'->new (">$f90")->print ($d->textContent);
    }
}

sub splitIntfb
{
  my $intfb =  shift;

  for my $openacc (0, 1)
    {
      my $suff = $openacc ? '_gpu' : '_cpu';
  
      my ($d) = &fxtran::parse (fragment => &slurp ($intfb));

      &addSuffix ($d, $suff, $openacc);

      (my $f = $intfb) =~ s/\.intfb.h$/$suff.intfb.h/;
      
      'FileHandle'->new (">$f")->print ($d->textContent);
    }
}

my ($file) = @ARGV;

if ($file =~ m/\.F90$/o)
  {
    &splitF90 ($file);
  }
elsif ($file =~ m/\.intfb\.h$/o)
  {
    &splitIntfb ($file);
  }



