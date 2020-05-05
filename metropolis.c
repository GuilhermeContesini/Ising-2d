#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define coup_A -1.0
#define coup_B 0.0
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void metropolis(long int size_int_, long int mcs_int_, int spin_ary_[size_int_], int dim_int_, int lat_ary_[dim_int_],double temp_dbl_, gsl_rng * rng_ptr_){
  int selct_spin_int=0;
  int trial_bool=1;
  int new_spin_int=0;
  int *ngbr_ptr=NULL;
  double dlt_ene_dbl=0.0; 
  double dlt_mag_dbl=0.0; 
  ngbr_ptr=malloc(2*dim_int_*sizeof(int));
  selct_spin_int=gsl_rng_uniform_int(rng_ptr_, size_int_);
  get_ngbr(selct_spin_int, dim_int_, ngbr_ptr, lat_ary_);
  new_spin_int=(-1)*spin_ary_[selct_spin_int];
  dlt_mag_dbl=(new_spin_int-spin_ary_[selct_spin_int]);
  dlt_ene_dbl=coup_A*dlt_mag_dbl*(spin_ary_[ngbr_ptr[0]]+spin_ary_[ngbr_ptr[1]]+spin_ary_[ngbr_ptr[2]]+spin_ary_[ngbr_ptr[3]])+coup_B*(dlt_mag_dbl);
  free(ngbr_ptr);
  if(dlt_ene_dbl>=0.0){
    trial_bool=gsl_rng_uniform(rng_ptr_)<exp(-dlt_ene_dbl/(temp_dbl_+1e-8))?1:0;
  }
  spin_ary_[selct_spin_int]=((trial_bool*new_spin_int))+ ((1-trial_bool)*spin_ary_[selct_spin_int]);
}