#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "gets.c"
#include "sets.c"
#include "metropolis.c"

int main(	int argc,	char *argv[]){
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// PROFILE
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	/*
	argv:
		sys_size_int_= INPUT
		meas_size_int= INPUT
		mcs_int_= INPUT
		lattice_size_x_int= INPUT
		lattice_size_y_int=sys_size_int_/lattice_size_x_int
		seed_int_=IPNUT
	*/
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// VARIABLES
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int i=0;
	int j=0;
	int sys_size_int=256;
	int meas_size_int=12;
	int dim_int=2;
	int lat_ary[dim_int];
	lat_ary[0]=16;
	lat_ary[1]=16;
	int * spin_ptr=NULL; 
	long int mcs_int=100000;
	long int mag_cor_time_int=0;
	long int ene_cor_time_int=0;
	long int seed_int=123456;
	// long int prob_distr_size_int=0;
	double init_temp_dbl=1.00;
	double end_temp_dbl=4.00;
	double dlt_temp_dbl=0.01;
	double temp_dbl=init_temp_dbl;
	// double *spin_prob_distr_ptr=NULL;
	double * meas_ene_ptr=NULL;
	double * meas_mag_ptr=NULL;
	double * meas_ene_sq_ptr=NULL;
	double * meas_mag_sq_ptr=NULL;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ARGV
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	sys_size_int=atoi(argv[1]);
	meas_size_int=atoi(argv[2]);
	mcs_int=atoi(argv[3]);
	lat_ary[0]=(int)atoi(argv[4]);
	lat_ary[1]=(int)(sys_size_int/lat_ary[0]);
	seed_int=atoi(argv[5]);
	printf("INPUT PARAMETERS:\n\tSystem size;\t\t%d\n\tMeasument size;\t\t%d\n\tMonte Carlo steps;\t%ld\n\tLattice dimensions;\t%d, %d\n\tRNG seed;\t\t%ld\n",sys_size_int,meas_size_int,mcs_int,lat_ary[0],lat_ary[1],seed_int);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// OUTPUT FILE
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	FILE *file_ptr=NULL;
	set_file(&file_ptr, sys_size_int, meas_size_int, mcs_int, seed_int);
	if(file_ptr==NULL){
		printf("FILE ERROR\n");
	}

	FILE *file_2_ptr=NULL;
	set_file_2(&file_2_ptr, sys_size_int, meas_size_int, mcs_int, seed_int);
	if(file_2_ptr==NULL){
		printf("FILE ERROR\n");
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// RNG
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	gsl_rng *rng_ptr=gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set( rng_ptr, seed_int);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MEM-ALLOCATION
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	spin_ptr=malloc(sys_size_int*sizeof(int));
	meas_ene_ptr=malloc(mcs_int*sizeof(double));
	meas_mag_ptr=malloc(mcs_int*sizeof(double));
	meas_ene_sq_ptr=malloc(mcs_int*sizeof(double));
	meas_mag_sq_ptr=malloc(mcs_int*sizeof(double));
	memset(spin_ptr, 0, sys_size_int);
	set_array_nil(sys_size_int, meas_ene_ptr);
	set_array_nil(sys_size_int, meas_mag_ptr);
	set_array_nil(sys_size_int, meas_ene_sq_ptr);
	set_array_nil(sys_size_int, meas_mag_sq_ptr);
	// prob_distr_size_int=pow(2, meas_size_int);
	// spin_prob_distr_ptr=malloc(prob_distr_size_int*sizeof(double));
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// INITIALIZATION
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	set_spin_rand(sys_size_int, spin_ptr, rng_ptr);
	// set_spin_des(sys_size_int, spin_ptr);
	// set_spin_coop(sys_size_int, spin_ptr);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// SIMULATION
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	printf("SIMUALTION STARTED...\n");
	while(temp_dbl<end_temp_dbl){
		for(i=0; i<mcs_int; i++){
			for(j=0; j<sys_size_int; j++){
				metropolis(sys_size_int, i, spin_ptr, dim_int, lat_ary, temp_dbl, rng_ptr);
			}
		  meas_ene_ptr[i]=get_ene(sys_size_int, spin_ptr, dim_int, lat_ary)/sys_size_int;
			meas_mag_ptr[i]=get_mag(sys_size_int, spin_ptr)/sys_size_int;
			meas_ene_sq_ptr[i]=get_ene_sq(sys_size_int, spin_ptr, dim_int, lat_ary)/sys_size_int;
			meas_mag_sq_ptr[i]=pow(get_mag(sys_size_int, spin_ptr)/sys_size_int,2);
			// printf("%.2lf\t%.2lf\t%.2lf\t%.2lf\n", meas_ene_ptr[i], meas_ene_sq_ptr[i], meas_mag_ptr[i], meas_mag_sq_ptr[i]);
		}
		// autocor_time_int=get_fft_cor_time(mcs_int, meas_ene_ptr);
		// mag_cor_time_int=get_int_cor_time(mcs_int, meas_mag_ptr);
		// ene_cor_time_int=get_int_cor_time(mcs_int, meas_ene_ptr);
		get_time_evo_samples(file_2_ptr, mcs_int, temp_dbl, meas_ene_ptr, meas_mag_ptr, meas_ene_sq_ptr, meas_mag_sq_ptr);
		get_samples(file_ptr, mcs_int, temp_dbl, ene_cor_time_int, mag_cor_time_int, meas_ene_ptr, meas_mag_ptr, meas_ene_sq_ptr, meas_mag_sq_ptr);
		temp_dbl=temp_dbl+dlt_temp_dbl;
	}
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MEMORY DEALLOCATION
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  free(spin_ptr);
	free(meas_ene_ptr);
	free(meas_mag_ptr);
	free(meas_ene_sq_ptr);
	free(meas_mag_sq_ptr);
  // free( spin_prob_distr_ptr ); 
	gsl_rng_free( rng_ptr );
  fclose ( file_ptr );
  fclose ( file_2_ptr );
  printf("SIMULATION FINISHED.\n");
  return 0;
}