#!/usr/bin/R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Library
library(ggplot2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }
	if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Correlation length FFT
cor_length_fft <- function( measure_points_array ){
	len_data <- length(measure_points_array)/2 
	half_measure_points_array <- measure_points_array[1:len_data]
	fft_measure_points_array <- fft(half_measure_points_array)
	normal_modes_float <- fft(fft_measure_points_array*Conj(fft_measure_points_array),inverse=TRUE)/len_data**2
	fft_cor_real <- (Re(normal_modes_float)-(mean(half_measure_points_array)**2) )/(var(half_measure_points_array )+1e-10 )
	return(fft_cor_real)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Specific Heat
specific_heat <- function( ene_ary, ene_sq_ary, temp_dbl ){
	data_len <- length( ene_ary )
	spcf_heat <- rep( 0, data_len )
	spcf_heat_stat <- rep( 0, 2 )
	for( i in seq( 1, data_len, 1 ) ) {
		spcf_heat[i] <- (1/(temp_dbl*temp_dbl+1e-8))*(ene_sq_ary[i]-(ene_ary[i]**2))/data_len
	}
	spcf_heat_stat[1] <- mean(spcf_heat)
	spcf_heat_stat[2] <- sd(spcf_heat)
	return( spcf_heat_stat )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Magnetic Susceptibility
mag_susc <- function( mag_ary, mag_sq_ary, temp_dbl ){
	data_len <- length( mag_ary )
	mag_qui <- rep( 0, data_len )
	mag_qui_stat <- rep( 0, 2 )
	for( i in seq( 1, data_len, 1 ) ) {
		mag_qui[i] <- (1/(temp_dbl+1e-8))*( mag_sq_ary[i]-(mag_ary[i]**2) )/data_len
	}
	mag_qui_stat[1] <- mean( mag_qui )
	mag_qui_stat[2] <- sd( mag_qui )
	return( mag_qui_stat )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Blocking_1
blocking_1 <- function( data_ary, ary_length ){
	division <- as.integer( ary_length/64 )
	ene_blocks <- split( data_ary, cut(seq_along(data_ary), division, labels=FALSE ))
	ene_blocks_avg_aux<-lapply(ene_blocks,mean)
	ene_blocks_avg<-array(as.numeric(unlist(ene_blocks_avg_aux)),dim=c(division))
	return( ene_blocks_avg )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Blocking_2
blocking_2 <- function( data_ary, ary_length ){
	division <- as.integer( 4*ary_length )
	ene_blocks <- split( data_ary, cut(seq_along(data_ary), division, labels=FALSE ))
	ene_blocks_avg <- array( as.numeric( unlist( ene_blocks_choosen_aux ) ), dim=c( division ) )
	ene_blocks_choosen_aux <- lapply( ene_blocks, sapply, sample(ene_blocks,1) )
	return( ene_blocks_avg )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Derivative
derivative <- function( data_ary, temp_ary, division ){
	data_ary_der <- rep( 0, length( data_ary )-1 )
	for( i in seq( 1, length( data_ary )-1, 1 ) ){
		data_ary_der[i] <- division*( data_ary[i+1] - data_ary[i] )/( temp_ary[i+1] - temp_ary[i] )
	}
	return( data_ary_der )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Error
error <- function(data_sd_ary){
	data_len <- length(data_sd_ary)
	return( qt(0.975,df=data_len-1)*data_sd_ary/sqrt(data_len) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Critical Temperature
critical_temp_1 <- function( data_ary ){
	
	# On this method I select the "max_samples" highest values and take the
	# average. In oder to accept this method, I make a assumption that the
	# function is symmetric and that the centrality is directly related to the
	# critical temperature.

	max_samples <- 10
	max_points_ary <- match( tail(data_ary[order(data_ary)], max_samples), data_ary ) 
	# cat( max_points_ary)
	return( mean(convert_index_temp[max_points_ary],trim=0.001) )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Critical Temperature
critical_temp_2 <- function( data_ary, temp_ary ){

	# On this method, firstly I invert the data so that the maximum value lies in
	# zero, in other words the maximum becomes the minimum. After we apply the
	# what I call the "Modefied" Newton method for finding the minimum, but
	# instead of calculating the value f(x) I calculate the first neighborhod
	# average, and instead of calculating the derivative f'(x) I use the first
	# neighborhod variance or standard deviation. The temperature that most
	# approximate zero will be the critical temperature.

	max_points <- max(data_ary)
	reverse_ary <- (max_points-data_ary)
	critical_temp <- 2.9 #Chute
	old_temperature <- 2.0 #Chute
	i <- which.min(abs(temp_ary-critical_temp))
	while(critical_temp != old_temperature){
		old_temperature <- critical_temp
		critical_temp <- temp_ary[i] - mean( c(reverse_ary[i+1],reverse_ary[i],reverse_ary[i-1]) )/sd(c(reverse_ary[i+1], reverse_ary[i], reverse_ary[i-1]))
		i <- which.min(abs(temp_ary-critical_temp))
	}
	return( critical_temp )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Critical Temperature
# critical_temp_3 <- function( data_ary ){

	# This method was inspired in the Fixed-point iteration method. Firstly the
	# method calculate the bigest and the lowest derivative, or the two biggest
	# derivative in absolute value. Here I assume or require the correct function
	# or distribution to be symmetric. Each derivative will give me a straight
	# line. The now is to find in which temperature the lines intercept, have the
	# same value. The temperature where these line intercept will be the critical
	# temperature.

# 	max_points <- max(data_ary)
# 	reverse_ary <- (max_points-data_ary)
# }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Critical Temperature
critical_temp_4 <- function(data_ary, temp_ary){
	# Once more, I invert the data so that the maximum value lies in zero, in
	# other words the maximum becomes the minimum. This is the secant method.
	# Instead os calculating the the analitical function every step I calculate
	# the nearst points on every interaction.

	max_points <- max(data_ary)
	reverse_ary <- (max_points-data_ary)

	critical_temp <- 2.9 #Chute
	old_temp <- 2.00 #Chute
	i <- which.min(abs(temp_ary-critical_temp))
	while(critical_temp != old_temp){
		old_temp <- critical_temp
		critical_temp <-temp_ary[i]-(temp_ary[i-1]-temp_ary[i])*( reverse_ary[i] )/( reverse_ary[i-1] - reverse_ary[i] )
		i <- which.min(abs(temp_ary-critical_temp))
	}
	return( critical_temp )
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Critical Temperature
# critical_temp_5 <- function(data_ary, temp_ary){

	# Teste de hipotese.

# }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main
ary_length <- 128
init_temp <- 4.00
end_temp <- 1.00
step_temp <- 0.01
division <- 100
temperature <- c(seq(init_temp*division,end_temp*division,-1)/division)
temp_len <- (init_temp-end_temp+step_temp)/step_temp

convert_index_temp <- rep(0,temp_len)
ene_avg_hist <- rep(0,temp_len)
ene_sd_hist <- rep(0,temp_len)
ene_sq_avg_hist <- rep(0,temp_len)
ene_sq_sd_hist <- rep(0,temp_len)
spcf_heat_avg <- rep(0,temp_len)
spcf_heat_sd <- rep(0,temp_len)
spcf_heat_list <- rep(0,2)

mag_avg_hist <- rep(0,temp_len)
mag_sd_hist <- rep(0,temp_len)
mag_sq_avg_hist <- rep(0,temp_len)
mag_sq_sd_hist <- rep(0,temp_len)
mag_qui_avg <- rep(0,temp_len)
mag_qui_sd <- rep(0,temp_len)
mag_qui_list <- rep(0,2)


temp <- init_temp
for(i in seq(1,temp_len)){
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	main_df <- read.table(
		paste("../Data/",ary_length*ary_length,"/time_series_",ary_length*ary_length,"-",i,".txt", sep=""),
		sep="\t",
		header=FALSE,
		skip="\n"
	)
	colnames(main_df) <- c("temperature", "mcs", "energy", "energy_sq", "magnetization", "magnetization_sq")
	convert_index_temp[i]<-temp
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	data_len <- length(main_df$energy)
	# ene_global_avg <- mean(main_df$energy)
	# ene_sq_global_avg <- mean(main_df$energy_sq)
	# mag_global_avg <- mean(main_df$magnetization)
	# mag_sq_global_avg <- mean(main_df$magnetization_sq)
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Blocking
	
	ene_blocks <- blocking_1( main_df$energy, ary_length )
	ene_avg_hist[i] <- mean( ene_blocks )
	ene_sd_hist[i] <- sd( ene_blocks )

	ene_sq_blocks <- blocking_1( main_df$energy_sq, ary_length )	
	ene_sq_avg_hist[i] <- mean( ene_sq_blocks )
	ene_sq_sd_hist <- sd( ene_sq_blocks )

	mag_blocks <- blocking_1( main_df$magnetization, ary_length )
	mag_avg_hist[i] <- mean( mag_blocks )
	mag_sd_hist[i] <- sd( mag_blocks )

	mag_sq_blocks <- blocking_1( main_df$magnetization_sq, ary_length )
	mag_sq_avg_hist[i] <- mean( mag_sq_blocks )
	mag_sq_sd_hist[i] <- sd( mag_sq_blocks )

	spcf_heat_list <- specific_heat( ene_blocks, ene_sq_blocks, temp)
	spcf_heat_avg[i] <- spcf_heat_list[1]
	spcf_heat_sd[i] <-  spcf_heat_list[2]

	mag_qui_list <- mag_susc( mag_blocks, mag_sq_blocks, temp)
	mag_qui_avg[i] <- mag_qui_list[1]
	mag_qui_sd[i] <- mag_qui_list[2]

	temp=temp-step_temp
	# cat(temp,"\n")
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Error
mag_avg_error <- error( mag_sd_hist )
ene_avg_error <- error( ene_sd_hist )
spcf_heat_error <- error( spcf_heat_sd )
mag_qui_error <- error( mag_qui_sd )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Critical Temperature Estimate

temp_1_near_crt_spcf_heat <- critical_temp_1( spcf_heat_avg )
temp_1_near_crt_mag_qui <- critical_temp_1( mag_qui_avg )

# temp_2_near_crt_spcf_heat <- critical_temp_2( spcf_heat_avg,temperature )
# temp_2_near_crt_mag_qui <- critical_temp_2( mag_qui_avg,temperature )

# temp_3_near_crt_spcf_heat <- critical_temp_3( spcf_heat_avg, temperature )
# temp_3_near_crt_mag_qui <- critical_temp_3( mag_qui_avg, temperature )

# temp_4_near_crt_spcf_heat <- critical_temp_4( spcf_heat_avg, temperature ) 
# temp_4_near_crt_mag_qui <-  critical_temp_4( mag_qui_avg, temperature ) 

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Critical Temperature via 1 method  mag suscpt:", temp_1_near_crt_mag_qui,'\n' )
cat("Critical Temperature via 1 method  spcf heat:", temp_1_near_crt_spcf_heat,'\n' )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cat("Critical Temperature via 2 method mag suscpt:", temp_2_near_crt_mag_qui,'\n' )
# cat("Critical Temperature via 2 method spcf heat:", temp_2_near_crt_spcf_heat,'\n' )

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cat("Critical Temperature via 3 method  mag suscpt:", temp_3_near_crt_mag_qui,'\n' )
# cat("Critical Temperature via 3 method  spcf heat:", temp_3_near_crt_spcf_heat,'\n' )

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cat("Critical Temperature via 4 method mag suscpt:", temp_4_near_crt_mag_qui,'\n' )
# cat("Critical Temperature via 4 method spcf heat:", temp_4_near_crt_spcf_heat,'\n' )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mag
mag_df <- data.frame(
	"temperature"=temperature,
	mag_avg_hist
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mag_plot <- ggplot(
	mag_df,
	aes(
		x=temperature,
		y=mag_avg_hist,
	)
)+ geom_point(
	aes(
		temperature,
		mag_avg_hist,
		color = "original_signal"
	),
	color='blue'
)+geom_errorbar(
	aes(
		ymin=mag_avg_hist-mag_avg_error,
		ymax=mag_avg_hist+mag_avg_error
	),
	color="black",
	width=0.05
)+ theme(
	text = element_text(size=10, family="LM Roman 10"),
	panel.border = element_rect(fill = NA),
	legend.box.background = element_rect(),
  legend.box.margin = margin(6, 6, 6, 6),
  panel.background = element_rect(fill = "snow3"),
  panel.grid.major = element_line(colour = "snow3"),
)+ ggtitle(
	"Magnetization x Temperature"
)+xlab("Temperature"
)+ylab("<M>"
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mag Qui
mag_qui_df <- data.frame(
	"temperature"=temperature,
	mag_qui_avg
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mag_qui_plot <- ggplot(
	mag_qui_df,
	aes(
		temperature,
		mag_qui_avg,
	)
)+geom_point(
	aes(
		temperature,
		mag_qui_avg,
		color = "original_signal"
	),
	color='red',
)+ggtitle(
	"Magnetic Susceptibility"
)+ theme(
	panel.border = element_rect(fill = NA),
	legend.box.background = element_rect(),
  legend.box.margin = margin(6, 6, 6, 6),
  panel.background = element_rect(fill = "snow3"),
  panel.grid.major = element_line(colour = "grey90"),
)+geom_errorbar(
	aes(
		ymin=mag_qui_avg-mag_qui_error,
		ymax=mag_qui_avg+mag_qui_error
	),
	color="black",
	width=0.05
)+ ggtitle(
	"Magnetic Susceptibility x Temperature"
)+xlab("Temperature"
)+ylab(c(expression(qui))
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mag Qui Derivative log-log
# mag_qui_der_df <- data.frame(
# 	"temperature"=temperature[-301],
# 	mag_qui_der
# )

# mag_qui_der_plot <- ggplot(
# 	mag_qui_der_df,
# 	aes(
# 		temperature,
# 		mag_qui_der,
# 	)
# )+geom_point(
# 	aes(
# 		temperature,
# 		mag_qui_der,
# 		color = "original_signal"
# 	),
# 	color='blue',
# )+ggtitle(
# 	"Magnetic Susceptibility"
# )+ theme(
# 	panel.border = element_rect(fill = NA),
# 	legend.box.background = element_rect(),
#   legend.box.margin = margin(6, 6, 6, 6),
#   panel.background = element_rect(fill = "snow3"),
#   panel.grid.major = element_line(colour = "grey90"),
# )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Energy
ene_df <- data.frame(
	"temperature"=temperature,
	ene_avg_hist
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ene_plot <- ggplot(
	ene_df,
	aes(
		temperature,
		ene_avg_hist,
	)
)+geom_point(
	aes(
		temperature,
		ene_avg_hist,
		color = "original_signal"
	),
	color='blue',
)+ggtitle(
	"Energy"
)+ theme(
	panel.border = element_rect(fill = NA),
	legend.box.background = element_rect(),
  legend.box.margin = margin(6, 6, 6, 6),
  panel.background = element_rect(fill = "snow3"),
  panel.grid.major = element_line(colour = "grey90"),
)+geom_errorbar(
	aes(
		ymin=ene_avg_hist-ene_avg_error,
		ymax=ene_avg_hist+ene_avg_error
	),
	color="black",
	width=0.05
)+ ggtitle(
	"Magnetic Susceptibility x Temperature"
)+xlab("Temperature"
)+ylab("<E>"
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Specific Heat
spcf_heat_df <- data.frame(
	"temperature"=temperature,
	spcf_heat_avg
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
spcf_heat_plot <- ggplot(
	spcf_heat_df,
	aes(
		temperature,
		spcf_heat_avg,
	)
)+geom_point(
	aes(
		temperature,
		spcf_heat_avg,
		color = "original_signal"
	),
	color='red',
)+ggtitle(
	"Specific Heat"
)+ theme(
	panel.border = element_rect(fill = NA),
	legend.box.background = element_rect(),
  legend.box.margin = margin(6, 6, 6, 6),
  panel.background = element_rect(fill = "snow3"),
  panel.grid.major = element_line(colour = "grey90"),
)+geom_errorbar(
	aes(
		ymin=spcf_heat_avg-spcf_heat_error,
		ymax=spcf_heat_avg+spcf_heat_error
	),
	color="black",
	width=0.05
)+ ggtitle(
	"Specific Heat"
)+xlab("Temperature"
)+ylab("c(T)"
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Specific Heat Derivative log-log
# spcf_heat_der_df <- data.frame(
# 	"temperature"=temperature[-301],
# 	spcf_heat_der
# )

# spcf_heat_der_plot <- ggplot(
# 	spcf_heat_der_df,
# 	aes(
# 		temperature,
# 		spcf_heat_der,
# 	)
# )+geom_point(
# 	aes(
# 		temperature,
# 		spcf_heat_der,
# 		color = "original_signal"
# 	),
# 	color='blue',
# )+ggtitle(
# 	"specific_heat Derivative in log-scale"
# )+ theme(
# 	panel.border = element_rect(fill = NA),
# 	legend.box.background = element_rect(),
#   legend.box.margin = margin(6, 6, 6, 6),
#   panel.background = element_rect(fill = "snow3"),
#   panel.grid.major = element_line(colour = "grey90"),
# )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggsave("magnetization.png",mag_plot)
ggsave("magnetic_susceptibility.png",mag_qui_plot)
ggsave("energy.png",ene_plot)
ggsave("specific_heat.png",spcf_heat_plot)

# setEPS()
# postscript("critical_phenomena.eps")


# multiplot(
# 	mag_plot,
# 	mag_qui_plot,
# 	# mag_qui_der_plot,
# 	ene_plot,
# 	spcf_heat_plot,
# 	# spcf_heat_der_plot,
# 	cols=2
# )

# dev.off()