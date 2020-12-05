## functions
# generating stimuli
IGen_fx <- function(m,V, feature_coding, firing_features){
	I <- lapply(V,function(i){rep(feature_coding[1],i)})
	# at each dimension i, turning on firing feature
	I <- lapply(c(1:m),function(i){
		stim_i <- I[[i]]
		stim_i[firing_features[i]] <- 1
		return(stim_i)
	})
	I <- unlist(I)
	return(I)
}

# equation (4)
eq4_fx <- function(I,i,H_j,V){
	end_pos <- sum(V[1:i])
	start_pos <- (end_pos-V[i])+1
	I_pos_i <- I[start_pos: end_pos]
	Hj_pos_i <- H_j[start_pos: end_pos]
	mue_i <- 0.5*sum(sqrt((I_pos_i-Hj_pos_i)^2))
	return(mue_i) # distance between cluster Hj and stimulus I at dimension i
}
mue_forEvery_DimAndCluster_fx <- function(j,H,I,m,V){
	H_j <- H[j,]
	mue_ij <- sapply(c(1:m),eq4_fx,I=I,H_j=H_j,V=V)
	return(mue_ij)
}
H_j_act_numerator_fx <- function(lambda_vec,params, mue_ij,j){
	numerator_j <- lambda_vec^params["r"]*exp(-lambda_vec* mue_ij[j,])
	return(numerator_j)
}
Eq4to7_fx <- function(H, M_HC, I_unlabeled, I, m, V, lambda_vec, maskedValues, params){
	
	# Equation 4
	if(is.null(H)){
	  H <- matrix(I, nrow=1) #CHANGE: Now uses I instead of I_unlabled (this is what Gureckis does)
	  new_recruit <- T
	}
	mue_ij <- t(sapply(c(1:nrow(H)),mue_forEvery_DimAndCluster_fx,H=H,I= I_unlabeled,m=m,V=V))
	
	# Equation 5
	H_j_act_numerator <- t(sapply(c(1:nrow(H)),H_j_act_numerator_fx, lambda_vec= lambda_vec,
	                              params=params, mue_ij= mue_ij))
	H_j_act_numerator <- rowSums(H_j_act_numerator,na.rm=T)
	#H_j_act_numerator <- rowSums(lambda_vec^params["r"]*exp(-lambda_vec* mue_ij),na.rm=T)
	H_j_act_denom <- sum(lambda_vec^params["r"],na.rm=T)
	H_j_act <- H_j_act_numerator/H_j_act_denom
	
	# Equation 6
	j_winner <- which(H_j_act==max(H_j_act)) # in case 2 or more clusters receive the same amount of activation
	j_winner <- j_winner[sample((1:length(j_winner)),1)]
	Act_after_competition <- (H_j_act^params["beta"]/sum(H_j_act^params["beta"]))*H_j_act
	H_j_out <- rep(0,nrow(H))
	H_j_out[j_winner] <- Act_after_competition[j_winner]
	
	# Equation 7
	# "When a cluster is recruited, weights from the unit to the output units are set to zero." (p. 316)
	
	#CHANGE: the previous matrix multiplication did not work
	if(is.null(M_HC)){M_HC <- matrix(rep(0,sum(V)),nrow=1)}
	C_zk_out <- M_HC[j_winner, maskedValues] * H_j_out[j_winner]

	return(list(H, mue_ij , M_HC, j_winner, H_j_out, C_zk_out))
}
# Equation 9
eq9_fx <- function(k, humbleTeacher, C_zk_out){
	if(humbleTeacher[k]==1){
		t_k <- 1
	}else{
		t_k <- 0
	}
	return(t_k)
}
# core function combining previous functions to the SUSTAIN-algorithm
trial_fx <- function( # fixed
						params, trial, block_stimuli, ResponseSet, maskedValues, maskedDim, m, V, 
						# for nested function Eq4to7, i.e., change recursively
							lambda_vec, H, M_HC,
								nStimCorrect, i
									){		
	I <- block_stimuli[trial,]
	I_unlabeled <- I
	I_unlabeled[maskedValues] <- 0
	
	Eq4to7 <- Eq4to7_fx( H=H, M_HC=M_HC, I_unlabeled=I_unlabeled, I=I, m=m, V=V, lambda_vec=lambda_vec,
	                     maskedValues= maskedValues, params= params)
	names(Eq4to7) <- c("H","mue_ij","M_HC","j_winner","H_j_out","C_zk_out")
	H <- Eq4to7$H
	M_HC <- Eq4to7$M_HC
	mue_ij <- Eq4to7$mue_ij
	j_winner <- Eq4to7$j_winner
	H_j_out <- Eq4to7$H_j_out
	C_zk_out <- Eq4to7$C_zk_out
	
	# Equation 8
	Pr_k <- exp(params["d"]* C_zk_out)/sum(exp(params["d"]* C_zk_out))
	Response <- sample(ResponseSet,1,prob= Pr_k)
	
	# Equation 9 (Forming target values t_zk)
	humbleTeacher <- I[maskedValues]
	t_zk <- sapply( c(1:V[maskedDim]), eq9_fx, humbleTeacher=humbleTeacher, C_zk_out=C_zk_out )
	correctResponse <- ResponseSet[which(humbleTeacher==1)]
	
	# Equation 10 (Deciding whether to recruit a new cluster or not)
	C_zk_out_max <- which(C_zk_out == max(C_zk_out))
	if(length(C_zk_out_max)>1){
	  t_value_at_Cmax <- 0
	  #print(C_zk_out)
	}else{
	  t_value_at_Cmax <- t_zk[C_zk_out_max]
	}
	
	if((i>1)&(t_value_at_Cmax!=1)){
		# "When a new cluster is recruited ..., it is centered on the misclassified input pattern... (p.316)"
		H <- rbind(H, I)
		rownames(H) <- paste("Cluster",c(1:nrow(H)))
		M_HC <- rbind(M_HC,rep(0,sum(V)))
		# "... and the clusters' activations and outputs are recalculated.
		#The new cluster then becomes the winner...(p.316)"
		Eq4to7 <- Eq4to7_fx( H=H, M_HC=M_HC, I_unlabeled=I_unlabeled, I=I, m=m, V=V, lambda_vec=lambda_vec,
		                     maskedValues= maskedValues, params=params)
		names(Eq4to7) <- c("H","mue_ij","M_HC","j_winner","H_j_out","C_zk_out")
		H <- Eq4to7$H
		M_HC <- Eq4to7$M_HC
		mue_ij <- Eq4to7$mue_ij
		j_winner <- Eq4to7$j_winner
		H_j_out <- Eq4to7$H_j_out
		C_zk_out <- Eq4to7$C_zk_out
		
		#CHANGE: Now also resets t values. This is significant for Eq #14 below
		humbleTeacher <- I[maskedValues]
		t_zk <- sapply( c(1:V[maskedDim]), eq9_fx, humbleTeacher=humbleTeacher, C_zk_out=C_zk_out )
	}

	# Equation 11 (applies only to unsupervised learning)

	# Equation 12 (Re-centering winning cluster)
	I_pos_ik <- I
	H_j_pos_ik <- H[j_winner,]
	delta_H_j_pos_ik <- params["eta"]*(I_pos_ik - H_j_pos_ik)
	H[j_winner,] <- H[j_winner,]+ delta_H_j_pos_ik

	# Equation 13 (Adjusting lambdas)
	eq_13_comp1 <- params["eta"]*exp(-lambda_vec*mue_ij[j_winner,])
	eq_13_comp2 <- 1-lambda_vec*Eq4to7$mue_ij[j_winner,]
	delta_lambda_i <- eq_13_comp1* eq_13_comp2
	lambda_vec <- c(lambda_vec+ delta_lambda_i)
	lambda_vec[length(lambda_vec)] <- 0
	#write(delta_lambda_i,"lambda.txt",append=T)
	
	# Equation 14 (Adjusting connections)
	delta_w_j_zk <- params["eta"]*(t_zk-C_zk_out)* H_j_out[j_winner]
	M_HC[j_winner, maskedValues] <- M_HC[j_winner, maskedValues]+ delta_w_j_zk
	#write(matrix(delta_w_j_zk,nrow=1),"delta_w_j_zk.txt",append=T)
	
	if(Response==correctResponse){nStimCorrect <- nStimCorrect+1}
		
	return(list(nStimCorrect, lambda_vec,H,M_HC))
	
}
block_fx <- function( params, agent, block, stimuli, ResponseSet, structure_name, 
							m, maskedDim, maskedValues, V, i,
								lambda_vec, H, M_HC
									 ){
		
		probOfError <- NULL # summarizes block performance
		
		block_stimuli <- stimuli [sample(nrow(stimuli)),]
		nStimCorrect <- 0
		
		for(trial in 1:nrow(block_stimuli)){
			
			i <- i+1
			
			trial_outcome <- trial_fx(params=params,trial=trial, block_stimuli= block_stimuli, ResponseSet = ResponseSet,
			                          maskedDim= maskedDim, maskedValues=maskedValues, m=m, V=V, lambda_vec= lambda_vec,
			                          H=H, M_HC= M_HC, nStimCorrect= nStimCorrect, i=i)
			names(trial_outcome) <- c("nStimCorrect","lambda_vec","H","M_HC")
			
			# for next trial
			lambda_vec <- trial_outcome$lambda_vec
			H <- trial_outcome$H
			M_HC <- trial_outcome$M_HC
			nStimCorrect <- trial_outcome$nStimCorrect
			
		}
		
		probOfError_block <- 1-(nStimCorrect/trial)
		
		block_performance_df <- data.frame( agent, structure_name, block, probOfError_block )
		
		block_list <- list(lambda_vec, H, M_HC, block_performance_df)
		names(block_list) <- c("lambda_vec", "H", "M_HC", "block_performance_df")
		
		return(block_list)
}


agent_fx <- function( params, agent, nBlocks, termCrit,
							labelQuery, ResponseSet, structure_name, stimulus_structure, 
								m, V, feature_coding
									){
	# "...32 blocks or ... 4 blocks without an error. ..., a block is defined as the presentation of each
  #item in a random order."(p. 319)
	
	stimuli <- t ( sapply ( stimulus_structure, IGen_fx, m=m, V=V, feature_coding=feature_coding ) )
	
	VLength <- sum(V)
	if(labelQuery==T){
		maskedDim <- m
		maskedValues <- (VLength-(V[m]-1)):VLength
	}else{
		maskedDim <- 1:(m-1)
		maskedValues <- 1:sum(V[1:(m-1)])
	}

	## INITIATING KNOWLEDGE STRUCTURE
	
	# "Initially, lambda_i is set to be broadly tuned with a value of 1.
	#The value of 1 is chosen because the maximal distance mue_ij is 1 and the optimal setting of lambda_i for
	#this case is 1 (i.e., Equation 13 equals zero). Under this scheme, lambda_i cannot become less than 1 but can
	#become more narrowly tuned." (p. 316)
	evolved_lambdas <- rep(1,m)
	# ??????? Does the model change its behavior if label-dimension is cancelled from calculating activation
	#values???????
	evolved_lambdas[maskedDim] <- NA
	# "Again SUSTAIN begins with a cluster centered on the first stimulus item." (p. 316)
	evolved_H <- NULL
	evolved_M_HC <- NULL
	
	errorFreeBlocks <- 0
	#agent_df <- NULL
	agent_performance_df <- NULL
	i <- 0
	
	for(block in 1: nBlocks){
		
		blockOutcome <- block_fx(params=params, agent=agent, block=block, stimuli=stimuli, ResponseSet=ResponseSet,
		                         structure_name= structure_name, m=m, maskedDim=maskedDim, maskedValues=maskedValues,
		                         lambda_vec=evolved_lambdas, H=evolved_H, M_HC=evolved_M_HC, i=i, V=V )
		
		agent_performance_df <- rbind(agent_performance_df, blockOutcome$block_performance_df)
		
		evolved_lambdas <- blockOutcome$lambda_vec
		evolved_H <- blockOutcome$H
		evolved_M_HC <- blockOutcome$M_HC
		
		if(blockOutcome$block_performance_df[,"probOfError_block"]==0){
		  errorFreeBlocks <- errorFreeBlocks+1
		}
		else{  #CHANGE: Paper speaks of 4 CONSECUTIVE error free blocks, not simply 4 error free blocks
		  errorFreeBlocks <- 0
		}
		if(errorFreeBlocks == termCrit){break}
	}
	
	return(agent_performance_df)
	
}
