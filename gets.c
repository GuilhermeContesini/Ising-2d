#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fft_complex.h>

#define coup_A -1.0
#define coup_B 0.0
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void get_ngbr(int i_, int dim_int_, int ngbr_ary_[2*dim_int_], int lat_ary_[dim_int_] ){
  int cur_x_int=i_%lat_ary_[ 0 ];
  int cur_y_int=i_/lat_ary_[ 0 ];
  int new_cur_x_int=0;
  int new_cur_y_int=0;
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Right
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  new_cur_x_int=(cur_x_int+1)%lat_ary_[0];
  new_cur_y_int=cur_y_int;
  ngbr_ary_[0]=new_cur_x_int+lat_ary_[0]*new_cur_y_int;
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Left
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  new_cur_x_int=(cur_x_int-1+lat_ary_[0])%lat_ary_[0];
  new_cur_y_int=cur_y_int;
  ngbr_ary_[1]=new_cur_x_int+lat_ary_[0]*new_cur_y_int;
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Up
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  new_cur_x_int=cur_x_int;
  new_cur_y_int=(cur_y_int-1+lat_ary_[1])%lat_ary_[ 1 ];
  ngbr_ary_[2]=new_cur_x_int+lat_ary_[0]*new_cur_y_int;
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Down 
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  new_cur_x_int=cur_x_int; 
  new_cur_y_int=(cur_y_int+1)%lat_ary_[1];
  ngbr_ary_[3]=new_cur_x_int+lat_ary_[0]*new_cur_y_int;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void get_label(int label_int_, int meas_size_int_, int spin_ary_[meas_size_int_]){
  int i=0;
  int aux_jnt=label_int_;

  for(i=0; i<meas_size_int_; i++){
    spin_ary_[i]=2*(aux_jnt%2)-1;
    aux_jnt=aux_jnt/2;
  }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_mean(int size_int_, double meas_ary_[size_int_]){
  int i=0;
  long double sum_dbl=0.0;
  for(i=0;i<size_int_;i++){
    sum_dbl=sum_dbl+ meas_ary_[i];
  }
  return sum_dbl/size_int_;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_mag(int size_int_, int spin_ary_[size_int_]){
  int i=0;
  double sum_dbl=0.0;
  for (i=0; i<size_int_; i++){
    sum_dbl=sum_dbl+spin_ary_[i];
  }
  return sum_dbl;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_ene(int size_int_, int spin_ary_[size_int_], int dim_int_, int lat_ary_[dim_int_]){
  int i=0;
  double ene_sum_dbl=0.0;
  int *ngbr_ptr=NULL;
  ngbr_ptr=malloc(2*dim_int_*sizeof(int));
  
  for(i=0; i<size_int_; i++){
    get_ngbr(i, dim_int_, ngbr_ptr, lat_ary_);
    ene_sum_dbl=ene_sum_dbl+(coup_A*(spin_ary_[i]*spin_ary_[ngbr_ptr[0]])+coup_A*(spin_ary_[i]*spin_ary_[ngbr_ptr[3]]))+coup_B*(spin_ary_[i]);
  }
  free(ngbr_ptr);
  return ene_sum_dbl;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_mag_sq(int size_int_, int spin_ary_[size_int_]){
  int i=0;
  double sum_dbl=0.0;
  for (i=0; i<size_int_; i++ ){
    sum_dbl=sum_dbl+(spin_ary_[i]*spin_ary_[i]);
  }
  return sum_dbl;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_ene_sq(int size_int_, int spin_ary_[size_int_], int dim_int_, int lat_ary_[dim_int_]){
  int i=0;
  double aux_ene_dbl=0.0;
  double ene_sum_dbl=0.0;
  int *ngbr_ptr=NULL;
  ngbr_ptr=malloc(2*dim_int_*sizeof(int));
  for(i=0; i<size_int_; i++){
    get_ngbr( i, dim_int_, ngbr_ptr, lat_ary_ );
    aux_ene_dbl=coup_A*(spin_ary_[i]*spin_ary_[ngbr_ptr[0]]+spin_ary_[i]*spin_ary_[ngbr_ptr[3]])+coup_B*(spin_ary_[i]);
    ene_sum_dbl=ene_sum_dbl+(aux_ene_dbl*aux_ene_dbl);
  }
  free(ngbr_ptr);
  return ene_sum_dbl;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_mag_qui(double temp_dbl_, double avg_mag_dbl_, double avg_mag_sq_dbl_ ){
  double mag_qui_dbl=0.0;
  mag_qui_dbl=(1/(temp_dbl_+1e-8))*(avg_mag_sq_dbl_-(avg_mag_dbl_*avg_mag_dbl_));
  return mag_qui_dbl;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double get_spcfc_heat(double temp_dbl_, double avg_ene_dbl_, double avg_ene_sq_dbl_){
  double spcfc_heat_dbl=0.0;
  spcfc_heat_dbl=(1/(temp_dbl_*temp_dbl_))*(avg_ene_sq_dbl_-(avg_ene_dbl_*avg_ene_dbl_));
  return spcfc_heat_dbl;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// double get_avg_ene(int meas_size_int_, int distr_size_int_, double prob_distr_int_[distr_size_int_] ){
//   int i=0;
//   int j=0;
//   int config_sum_int=pow(2, meas_size_int_);
//   double config_ene_dbl=0.0;
//   double avg_ene_dbl=0.0;

//   for(j=0;j<config_sum_int; j++){
//     config_ene_dbl=0.0;
//     int *config_ary=NULL;
//     config_ary=malloc(meas_size_int_*sizeof(int));
//     get_label(j, meas_size_int_, config_ary);

//     for(i=0; i<meas_size_int_-1; i++){
//       config_ene_dbl=config_ene_dbl+coup_A*(config_ary[i]*config_ary[i+1])+coup_B*(config_ary[i]);
//     }
//     config_ene_dbl=config_ene_dbl+coup_A*(config_ary[i+1]*config_ary[0])+coup_B*(config_ary[i+1]);
//     avg_ene_dbl=avg_ene_dbl+(config_ene_dbl*prob_distr_int_[j]);
//     free(config_ary);
//   }
//   return avg_ene_dbl/meas_size_int_;
// }
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// double get_avg_entropy(int distr_size_int_, double prob_distr_int_[ distr_size_int_ ] ){
//   int i=0;
//   double s_dbl=0.0;

//   for(i=0; i<distr_size_int_; i++){
//     s_dbl=s_dbl-(prob_distr_int_[i]*log(1.e-6+prob_distr_int_[i]) );
//   }
//   return s_dbl;
// }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// int get_fft_cor_time(int size_int_, double meas_ary_[size_int_]){
//   int i=0;
//   int cor_time_int=0;
//   gsl_complex_packed_array orig_ary=meas_ary_;
//   gsl_complex_packed_array conj_ary=;
//   gsl_complex_packed_array fft_ary_ary=meas_ary_;
//   size_t stride_int=1;
//   size_t n_int=size_int_;
//   double complex *fft_meas_ary=NULL;
//   fft_meas_ary=malloc(size_int_*sizeof(double complex));
//   memset(fft_meas_ary, 0, size_int_);
//   gsl_fft_complex_radix2_forward(orig_ary, stride_int, n_int);
//   conj_ary=conj(orig_ary);
//   for(i = 0; i < size_int_; ++i)
//   {
//     fft_ary_ary[i]=(orig_ary[i]*conj_ary[i])/128;
//   }
//   gsl_fft_complex_radix2_inverse(fft_ary_ary,stride_int,n_int);


//   return cor_time_int;
// }
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// int get_clas_cor_time(int size_int_, double meas_ary_[size_int_]){
//   int i=0;
//   int j=0;
//   int cor_length_int=0;
//   double mean_dbl=0.0;
//   double *sum_ary=NULL;
//   sum_ary=malloc(size_int_*sizeof(double));
//   memset(sum_ary, 0, size_int_/2);
//   mean_dbl=get_mean(size_int_,meas_ary_);
//   for(j=0; j<size_int_/2;j++){
//     for(i=0; i<size_int_/2; ++i){
//       sum_ary[i]=sum_ary[i]+(meas_ary_[i]*meas_ary_[i+j]);
//     }
//   }
//   free(sum_ary);
//   return cor_length_int;
// }
// cor(a[])=1/n 1/((a**2)-(a)**2) sum i=0 n/2 sum j=0 n/2 a[i]*a[i+j] -(a[])**2
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void get_time_evo_samples(FILE *file_ptr, int mcs_int_, double temp_dbl_, double meas_ene_ary_[mcs_int_], double meas_mag_ary_[mcs_int_], double meas_ene_sq_ary_[mcs_int_], double meas_mag_sq_ary_[mcs_int_]){
  int i=0;
  for( i=0; i<mcs_int_; i++) {
    // printf("%.2lf\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",temp_dbl_, i, meas_ene_ary_[i], meas_ene_sq_ary_[i], meas_mag_ary_[i], meas_ene_sq_ary_[i]);
    fprintf(file_ptr, "%.2lf\t%d\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",temp_dbl_, i, meas_ene_ary_[i], meas_ene_sq_ary_[i], meas_mag_ary_[i], meas_mag_sq_ary_[i]);
  }
  fprintf(file_ptr,"\n");
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void get_samples(FILE *file_ptr, int mcs_int_, double temp_dbl_, double ene_cor_time_int_, double mag_cor_time_int_, double meas_ene_ary_[mcs_int_], double meas_mag_ary_[mcs_int_], double meas_ene_sq_ary_[mcs_int_], double meas_mag_sq_ary_[mcs_int_]){
  double avg_ene_dbl=0.0;
  double avg_mag_dbl=0.0;
  double avg_ene_sq_dbl=0.0;
  double avg_mag_sq_dbl=0.0;
  double spcfc_heat_dbl=0.0;
  double mag_qui_dbl=0.0;
  avg_ene_dbl=get_mean(mcs_int_,meas_ene_ary_);
  avg_mag_dbl=get_mean(mcs_int_,meas_mag_ary_);
  avg_ene_sq_dbl=get_mean(mcs_int_,meas_ene_sq_ary_);
  avg_mag_sq_dbl=get_mean(mcs_int_,meas_mag_sq_ary_);
  spcfc_heat_dbl=get_spcfc_heat(temp_dbl_, avg_ene_dbl, avg_ene_sq_dbl);
  mag_qui_dbl=get_mag_qui(temp_dbl_, avg_mag_dbl, avg_mag_sq_dbl);
  // printf("%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",temp_dbl_, avg_ene_dbl, avg_ene_sq_dbl, spcfc_heat_dbl, avg_mag_dbl, avg_mag_sq_dbl, mag_qui_dbl);
  fprintf(file_ptr, "%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", temp_dbl_, avg_ene_dbl, avg_ene_sq_dbl, spcfc_heat_dbl, avg_mag_dbl, avg_mag_sq_dbl, mag_qui_dbl);
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
