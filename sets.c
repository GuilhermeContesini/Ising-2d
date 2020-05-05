#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_file(FILE **file_ptr_, int sys_size_int, int meas_size_int, int mcs_int, int seed_int_){
	char *begin_str="./Data/data_file_Size-"; 
	char *btw_str="_Meas-";
	char *mcs_str="_MCS-";
	char *seed_str="_Seed-";
	char *ext_str=".txt";
	char file_str[strlen(begin_str)+strlen(ext_str)+strlen(btw_str)+strlen(seed_str)+strlen(mcs_str)+20];
	snprintf(file_str, sizeof(file_str), "%s%d%s%d%s%d%s%d%s", begin_str, sys_size_int, btw_str, meas_size_int, mcs_str, mcs_int, seed_str, seed_int_, ext_str );
	*file_ptr_=fopen(file_str, "w");
	if(*file_ptr_==NULL){
		printf("FILE ERROR!");      
	}
	printf("DATA FILE OPEN!\n");
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_file_2(FILE **file_ptr_, int sys_size_int, int meas_size_int, int mcs_int, int seed_int_){
	char *begin_str="./Data/data_time_series_Size-"; 
	char *btw_str="_Meas-";
	char *mcs_str="_MCS-";
	char *seed_str="_Seed-";
	char *ext_str=".txt";
	char file_str[strlen(begin_str)+strlen(ext_str)+strlen(btw_str)+strlen(seed_str)+strlen(mcs_str)+20];
	snprintf(file_str, sizeof(file_str), "%s%d%s%d%s%d%s%d%s", begin_str, sys_size_int, btw_str, meas_size_int, mcs_str, mcs_int, seed_str, seed_int_, ext_str );
	*file_ptr_=fopen(file_str, "w");
	if(*file_ptr_==NULL){
		printf("FILE ERROR!");      
	}
	printf("TIME SERIES FILE OPEN!\n");
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_spin_rand(int size_int_, int spin_ary_[size_int_], gsl_rng * rng_gsl_ptr_){
	int i=0;
	for( i=0; i<size_int_; i++ ){
		spin_ary_[i]=-1;
		if( gsl_rng_uniform( rng_gsl_ptr_ ) < 0.5 ){
			spin_ary_[i]=1;
		}	
	}
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_spin_des(int size_int_, int spin_ary_[size_int_]){
  int i=0;
  for(i=0; i < size_int_; i++){
    spin_ary_[i]=-1;
  }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_spin_coop(int size_int_,int spin_ary_[size_int_]){
  int i=0;
  for(i=0; i<size_int_; i++){
		spin_ary_[i]=1;
  }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_array_nil(int size_int_,double spin_ary_[size_int_]){
  int i=0;
  for(i=0; i<size_int_; i++){
		spin_ary_[i]=0.0;
  }
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void set_prob_distr( int size_int_, int spin_ary_[ size_int_ ], int meas_size_int_, int prob_distr_size_int_, double prob_distr_int_[prob_distr_size_int_]){
	int ary_i=0;
	int j=0;
	int label_int=0;
	double sum_double=0.0;
  double sum_prob_ary[ prob_distr_size_int_ ];
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Array-Counter Cleaner 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(ary_i=0; ary_i<=prob_distr_size_int_; ary_i++){
		sum_prob_ary[ary_i]=0.0;
		prob_distr_int_[ary_i]=0.0;
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Count which and how many micro states there are.
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(ary_i=0; ary_i<=size_int_-meas_size_int_; ary_i++){
		label_int=0;
		for(j=meas_size_int_-1; j>=0; j--){
			label_int=+(((1+spin_ary_[ary_i+j])/2)*pow( 2, j ));
		}
		sum_prob_ary[ label_int ]++;
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Normalization of Probability Distribution .
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(ary_i=0; ary_i<=prob_distr_size_int_; ary_i++ ){
		sum_double=+sum_prob_ary[ary_i];
	}
	for(ary_i=0; ary_i < prob_distr_size_int_; ary_i++	){
		prob_distr_int_[ary_i]=sum_prob_ary[ary_i]/sum_double; 
	}
}
// ainda falta a condição de contorno !