

Vfamilies<-function(dat, depth, NP, N){
	
	temp<-matrix(0, ncol=3, nrow=depth)
	for(z in 1:depth){
		a=0
		d=0
		if(z!=0){a=1}
		if(a==1){d=z-1}
		
		first=a*np+d*N + 1
		last= np + z*N
		Dat=as.data.frame(dat[first:last,])
		Dat[,2]<-as.factor(Dat[,2])
		Dat[,3]<-as.factor(Dat[,3])
		Table <- summary(aov(Gm ~ 1+SIRE+(SIRE:DAM),data=Dat))
		DF <- as.matrix(Table[[1]][1]$"Df")
		SS <- as.matrix(Table[[1]][2]$"Sum Sq")
		n_sire <- length(levels(Dat[,3]))
		n_dam <- length(levels(Dat[,2]))/n_sire
		n_obs <- nrow(Dat)/(n_sire*n_dam)
		Ve <- SS[3,]/DF[3,]
		Vd <- (SS[2,1]/DF[2,1]-Ve)/n_obs    # dam var comp
		Vs <- (SS[1,1]/DF[1,1]-Ve-n_obs*Vd)/(n_dam*n_obs)   # sire var comp
		temp[z,]<-c(Ve,Vd,Vs)
		
		
	}
		
	return(temp)
}

StrassenCprod <-function(A, B){
	dimA <- dim(A)
	dimB <- dim(B)
	if(is.null(dimA)==TRUE | is.null(dimB)==TRUE){stop}
	if(dimA[1]==dimA[2] & dimA[1]==dimB[1] & dimB[1]==dimB[2]){
		a <- dimA[1]/2
		A11 <- A[1:a,1:a]			;	A12 <- A[1:a,(a+1):(2*a)]
		A21 <- A[(a+1):(2*a), 1:a]  ;   A22 <- A[(a+1):(2*a),(a+1):(2*a)]

		B11 <- B[1:a,1:a]			;	B12 <- B[1:a,(a+1):(2*a)]
		B21 <- B[(a+1):(2*a), 1:a]  ;   B22 <- B[(a+1):(2*a),(a+1):(2*a)]

		M1 <- crossprod((A11 + A22), (B11 + B22))
		M2 <- crossprod((A21 + A22), B11)
		M3 <- crossprod(A11, (B12 - B22))
		M4 <- crossprod(A22, (B21 - B11))
		M5 <- crossprod((A11 + A22), B22)
		M6 <- crossprod((A21 - A11), (B11 + B12))
		M7 <- crossprod((A12 - A22), (B21 + B22))
		
		C11 <- M1 + M4 - M5 + M7
		C12 <- M3 + M5
		C21 <- M2 + M4
		C22 <- M1 - M2 + M3 + M6
		
		C <- cbind(rbind(C11, C12), rbind(C21, C22))
		return(C)
		
	}

}



# ================================================================================================
# function to set phenotypic (co)variance matrices for traits and their joint covariance section
# ind, FS mean, HS mean, dam I, sire I, mean of half-sibs' dams I
P_matr <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,6,6)

	CovamD <- 4*Covamd # dam cov from dam family cov
	Vp <- Vas+Vad+Va0/2+VmD+CovamD+Vc+Ve; CovFp <- Vas + Vad + VmD + CovamD
	VFS <- Vas+Vad+VmD+CovamD+Vc+(Va0/2+Ve)/(n-1) 
	CovHp <- Vas; VHS <- Vas+(Vad+VmD+CovamD+Vc)/d +(Va0/2+Ve)/(d*n)
	CovId <- (1-kd)*(CovaI/2+CovmI); CovIs <- (1-ks)*(CovaI/2)	

	#   ind            FS mean          HS mean           dam index         sire index     dam index mean
	PP[1,1]<-Vp;     PP[1,2]<-CovFp;  PP[1,3]<-CovHp;  PP[1,4]<-CovId;    PP[1,5]<-CovIs;    PP[1,6]<- 0  
	PP[2,1]<-CovFp;  PP[2,2]<-VFS;    PP[2,3]<-CovHp;  PP[2,4]<-CovId;    PP[2,5]<-CovIs;    PP[2,6]<- 0  
	PP[3,1]<-CovHp;  PP[3,2]<-CovHp;  PP[3,3]<-VHS;    PP[3,4]<- 0;       PP[3,5]<-CovIs;    PP[3,6]<- CovId/d
	PP[4,1]<-PP[1,4];PP[4,2]<-PP[2,4];PP[4,3]<-PP[3,4];PP[4,4]<-(1-kd)*VI;PP[4,5]<-0;        PP[4,6]<- 0
	PP[5,1]<-PP[1,5];PP[5,2]<-PP[2,5];PP[5,3]<-PP[3,5];PP[5,4]<- 0;       PP[5,5]<-(1-ks)*VI;PP[5,6]<- 0 
	PP[6,1]<-PP[1,6];PP[6,2]<-PP[2,6];PP[6,3]<-PP[3,6];PP[6,4]<- 0;       PP[6,5]<-0;        PP[6,6]<-(1-kd)*VI/d 
	return(PP)
	
}

P_matr_a <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,6,6)

	CovamD <- 4*Covamd # dam cov from dam family cov
	Vp <- Vas+Vad+Va0/2+VmD+CovamD+Vc+Ve; CovFp <- Vas + Vad + VmD + CovamD +Vc
	VFS <- Vas+Vad+VmD+CovamD+Vc+(Va0/2+Ve)/(n-1) 
	CovHp <- Vas; VHS <- Vas+(Vad+VmD+CovamD+Vc)/d +(Va0/2+Ve)/(d*n)
	CovId <- (1-kd)*(CovaI/2); CovIs <- (1-ks)*(CovaI/2)	

	#   ind            FS mean          HS mean           dam index         sire index     dam index mean
	PP[1,1]<-Vp;     PP[1,2]<-CovFp;  PP[1,3]<-CovHp;  PP[1,4]<-CovId;    PP[1,5]<-CovIs;    PP[1,6]<- 0  
	PP[2,1]<-CovFp;  PP[2,2]<-VFS;    PP[2,3]<-CovHp;  PP[2,4]<-CovId;    PP[2,5]<-CovIs;    PP[2,6]<- 0  
	PP[3,1]<-CovHp;  PP[3,2]<-CovHp;  PP[3,3]<-VHS;    PP[3,4]<- 0;       PP[3,5]<-CovIs;    PP[3,6]<- CovId/d
	PP[4,1]<-PP[1,4];PP[4,2]<-PP[2,4];PP[4,3]<-PP[3,4];PP[4,4]<-(1-kd)*VI;PP[4,5]<-0;        PP[4,6]<- 0
	PP[5,1]<-PP[1,5];PP[5,2]<-PP[2,5];PP[5,3]<-PP[3,5];PP[5,4]<- 0;       PP[5,5]<-(1-ks)*VI;PP[5,6]<- 0 
	PP[6,1]<-PP[1,6];PP[6,2]<-PP[2,6];PP[6,3]<-PP[3,6];PP[6,4]<- 0;       PP[6,5]<-0;        PP[6,6]<-(1-kd)*VI/d 
	return(PP)
	
}

P_matr_index <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,3,3)

	CovamD <- 4*Covamd # dam cov from dam family cov
	Vp <- Vas+Vad+Va0/2+VmD+CovamD+Vc+Ve; CovFp <- Vas + Vad + VmD + CovamD 
	VFS <- Vas+Vad+VmD+CovamD+Vc+(Va0/2+Ve)/(n-1) 
	CovHp <- Vas; VHS <- Vas+(Vad+VmD+CovamD+Vc)/d +(Va0/2+Ve)/(d*n)
		

	#   ind            FS mean          HS mean     
	PP[1,1]<-Vp;     PP[1,2]<-CovFp;  PP[1,3]<-CovHp   
	PP[2,1]<-CovFp;  PP[2,2]<-VFS;    PP[2,3]<-CovHp    
	PP[3,1]<-CovHp;  PP[3,2]<-CovHp;  PP[3,3]<-VHS    
	return(PP)
	
}


P_matr_ph <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks){

	PP <- matrix(0,1,1)
	
	CovamD <- 4*Covamd
	Vp <- Vas+Vad+Va0/2+VmD+CovamD+Vc+Ve
	
	
	PP[1,1]<-Vp
	return(PP)
}		


P_matr_independent <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,3,3)
	
	CovamD <- 4*Covamd # dam cov from dam family cov
	
	VpMinusFS <- ((n-1)*(Va0/2+Ve))/(n)
	VFSminusHS <- (d-1)*(Vad + VmD + CovamD)/(d) + (d*n-(n-1))*(Va0/2+Ve)/(d*n*(n-1))
	VHS <- Vas + (Vad + VmD + CovamD)/(d) +(Va0/2+Ve)/(d*n)
	

	#   ind            FS mean          HS mean     
	PP[1,1]<-VpMinusFS;     	PP[1,2]<-0;  					PP[1,3]<-0   
	PP[2,1]<-0; 		 		PP[2,2]<-VFSminusHS;    		PP[2,3]<-0    
	PP[3,1]<-0;			  		PP[3,2]<-0;  					PP[3,3]<-VHS    
	return(PP)
	
}


P_matr_invMatingR <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,3,3)
	
	CovamD <- 4*Covamd # dam cov from dam family cov
	
	VpMinusFS <- ((n-1)*(Va0/2+Ve))/(n)
	VFSminusHS <- ((d-1)/d)*(Vas + ((Va0/2+Ve)/(n)))
	VHS <- Vad + VmD + CovamD + Vas/((d-1))+(1/(n*(d-1)))*(Va0/2+Ve)
	

	#   ind            FS mean          HS mean     
	PP[1,1]<-VpMinusFS;     	PP[1,2]<-0;  					PP[1,3]<-0   
	PP[2,1]<-0; 		 		PP[2,2]<-VFSminusHS;    		PP[2,3]<-0    
	PP[3,1]<-0;			  		PP[3,2]<-0;  					PP[3,3]<-VHS    
	return(PP)
	
}


P_matr_independent_reindeer <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) { 
  
	PP <- matrix(0,2,2)
	
	CovamD <- 4*Covamd # dam cov from dam family cov
	
	VpMinusHS <- (d-1)*(Vas + Va0/2 + Ve)/(d)
	
	VHS <- Vad + VmD + CovamD + (Va0/2+Ve+Vas)/d 
	

	#   ind            FS mean          HS mean     
	PP[1,1]<-VpMinusHS;     	PP[1,2]<-0  					  
	PP[2,1]<-0; 		 		PP[2,2]<-VHS    		   
	    
	return(PP)
	
	
}


P_matr_independent_Wparents <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,6,6)
	CovamD <- 4*Covamd # dam cov from dam family cov
	
	VpMinusFS <- ((n-1)*(Va0/2+Ve))/n
	VFSminusHS <- ((d-1)*(Vad + VmD + CovamD + ((Va0/2+Ve)/n))/d)
	VHS <- Vas + (Vad + VmD + CovamD)/d +(Va0/2+Ve)/(d*n)
	CovId <- (1-kd)*(CovaI/2+CovmI); CovIs <- (1-ks)*(CovaI/2)	

	#   ind            FS mean          HS mean           dam index         sire index     dam index mean
	PP[1,1]<-VpMinusFS;     PP[1,2]<-0;  			PP[1,3]<-0;  		PP[1,4]<- 0;    	PP[1,5]<-0;   		PP[1,6]<- 0  
	PP[2,1]<-0;  			PP[2,2]<-VFSminusHS;   	PP[2,3]<-0;  		PP[2,4]<-CovId;    	PP[2,5]<-0;    		PP[2,6]<- 0  
	PP[3,1]<-0;  			PP[3,2]<-0;  			PP[3,3]<-VHS;   	PP[3,4]<- 0;       	PP[3,5]<-CovIs;    	PP[3,6]<- CovId/d
	PP[4,1]<-PP[1,4];		PP[4,2]<-PP[2,4];		PP[4,3]<-PP[3,4];	PP[4,4]<-(1-kd)*VI;	PP[4,5]<-0;        	PP[4,6]<- 0
	PP[5,1]<-PP[1,5];		PP[5,2]<-PP[2,5];		PP[5,3]<-PP[3,5];	PP[5,4]<- 0;       	PP[5,5]<-(1-ks)*VI;	PP[5,6]<- 0 
	PP[6,1]<-PP[1,6];		PP[6,2]<-PP[2,6];		PP[6,3]<-PP[3,6];	PP[6,4]<- 0;       	PP[6,5]<-0;        	PP[6,6]<-(1-kd)*VI/d 
	return(PP)
	
}




P_matr_Grandsire <- function(d,n,Va0,Vas,Vad,Ve,VmD,Covamd, Covam0, Vc,CovaI,CovmI,VI,kd,ks) {
  
	PP <- matrix(0,2,2)
	TT<-1
	CovamD <- 4*Covamd # dam cov from dam family cov
	VAS<-4*Vas
	VAD<-4*Vad
	d2<-0.5*d*n
	
	Progeny <- VAS/4 + (VAD/4 + VmD + CovamD)/d + (Va0/2 + Ve)/(d*n)
	GrandProgeny <- VAS/16 + VmD/4 + CovamD/4 + ((3/16)*VAD + (3/4)*VmD + (3/4)* CovamD)/d + (VAD/4)/(d2) + (Va0/2 + Ve)/(d2*n)
	
	CovPGp <- VAS/8  + CovamD/4 + (VAD/8 + CovamD/4 + VmD/2)/(d) + (Va0/4 + CovamD/2)/(d*n*0.5)
	#     
	PP[1,1]<-Progeny;     		PP[1,2]<-CovPGp  					  
	PP[2,1]<-CovPGp; 		 	PP[2,2]<-GrandProgeny    		    
	   
	return(PP)
	
}



					
G_matr <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	g_a <- c(Vas+Vad+Va0/2+Covam_D/2,				Vas+Vad+Covam_D/2,					Vas,						(1-kd)*CovaI/2,		(1-ks)*CovaI/2, 0 )
	g_m <- c(Covams+Covamd+((1/2)*Covam0)+VmD/2, 	Covams+Covamd+VmD/2, 				Covams,						(1-kd)*CovmI/2,		(1-ks)*CovmI/2, 0 )       # cov(m,HS)=0 as ind's dam and half-sibs' dams unrelated
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}

G_matr_a <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	g_a <- c(Vas+Vad+Va0/2,										Vas+Vad,					Vas,					(1-kd)*CovaI/2,		(1-ks)*CovaI/2, 0 )
	g_m <- c(Covams+Covamd+((1/2)*Covam0)+VmD/2, 	Covams+Covamd+VmD/2, 					Covams,					(1-kd)*CovmI/2,		(1-ks)*CovmI/2, 0 )       # cov(m,HS)=0 as ind's dam and half-sibs' dams unrelated
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}

G_matr_index <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	g_a <- c(Vas+Vad+Va0/2+Covam_D/2,				Vas+Vad+Covam_D/2,					Vas )
	g_m <- c(Covams+Covamd+((1/2)*Covam0)+VmD/2, 	Covams+Covamd+VmD/2, 				Covams)    
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}

#test for selecting direct effect only
G_matr_index_a <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	g_a <- c(Vas+Vad+Va0/2,							Vas+Vad,							Vas)
	g_m <- c(Covams+Covamd+((1/2)*Covam0)+VmD/2, 	Covams+Covamd+VmD/2, 				Covams)    
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}



G_matr_ph <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	g_a <- c(Vas+Vad+Va0/2+Covam_D/2)
	g_m <- c(Covams+Covamd+((1/2)*Covam0)+VmD/2)       # cov(m,HS)=0 as ind's dam and half-sibs' dams unrelated
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}



G_matr_independent <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	
	g_a <- c(Va0/2,				Vad+Covam_D/2,					Vas)
	g_m <- c(Covam0/(2), 		(Covam_D/4+VmD/2), 				Covams)  
	
	GG <- cbind(g_a, g_m)
 
	return(GG)
}

G_matr_invMatingR <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	
	g_a <- c(Va0/2,				Vas,					Vad+Covam_D/2)
	g_m <- c(Covam_D/(2), 		Covam_D/4, 				Covam_D/4 + VmD/2)    
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}



G_matr_independent_reindeer <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	Covam_S <- 4*Covams
	g_a <- c(Vas + Va0/2, 		Vad+Covam_D/2)
	g_m <- c(Covams,			Covam_D/4 + VmD/2)  
	
	GG <- cbind(g_a, g_m)
  
	return(GG)
}



G_matr_independent_Wparents <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	g_a <- c((n-1)*(Va0/2)/n,				(d-1)*(Vad + Covam_D/2 + Va0/(2*n))/d,						Vas + (Vad + Covam_D/2)/d + Va0/(2*d*n),						(1-kd)*CovaI/2,		(1-ks)*CovaI/2, 0 )
	g_m <- c((n-1)*Covam_D/(2*n), 			(d-1)*(Covam_D/4+ VmD/2 + Covam_D/(2*n))/d, 				Covam_D/4 + (Covam_D/4 + VmD/2)/d + Covam_D/(2*d*n),			(1-kd)*CovmI/2,		(1-ks)*CovmI/2, 0 )       # cov(m,HS)=0 as ind's dam and half-sibs' dams unrelated
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}

#inversed mating ratio
G_matr_invMatingR <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	
	g_a <- c(Va0/2,				Vas,					Vad+Covam_D/2)
	g_m <- c(Covam_D/(2), 		Covam_D/4, 				Covam_D/4 + VmD/2)    
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}


G_matr_Grandsire <- function(d,n,Va0,Vas,Vad,Vm0,VmD,Covam0,Covams,Covamd,CovaI,CovmI,kd,ks,VI) {
 
	Covam_D <- 4*Covamd # covar among dams (computed from dam family covariance)
	VAS<-4*Vas
	VAD<-4*Vad
	
	g_a <- c(VAS/2,				VAS/4+Covam_D/2)
	g_m <- c(Covam_D/(2), 		Covam_D/4 + VmD/2)    
 
	GG <- cbind(g_a, g_m)
  
	return(GG)
}



########################################Function to calculate ################################# 
	calK<-function(p){
		x <- qnorm(p,lower.tail=F)
		z <- dnorm(x)
		i <- z/p  # selection intensity
		k <- i*(i-x)
		ik<-c(i,k)
		return(ik)
	}

#########################################################################################


PseudoBLUP<-function(ND=0, NFS=0, Gen=1, Va=0, Vm=0, Covam=0, Vc=0, Ve=0, w=c(1,1), ps=0, pd=0, iter=10, 
PhenoSelect=FALSE, IndexSelect=FALSE, Inbr=0, DirectOnly=FALSE, ws=c(1,1), independent=FALSE, independent_reindeer=FALSE, inversedMatingRatio=FALSE, BLUP=FALSE,
InDbSel=c(TRUE, TRUE, TRUE), independent2=FALSE, InDbSel2=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), InDbRein=c(TRUE,TRUE), Grandsire=FALSE){
########################################Phenotypic selection###############################
	
	ram <- Covam/sqrt(Va*Vm)
	DirectV=NULL
	DirectResponse=NULL
	DirectRaa=NULL
	MaternalV=NULL
	MaternalResponse=NULL
	MaternalRmm=NULL
	resultsSUM=NULL
	Covariance=NULL
	DirectRaI=NULL
	MaternalRmI=NULL
	IndexVar=NULL
	RAM=NULL
	Qresp=c(0,0)
	#Pvar=
	
	w<-ws
	if(DirectOnly==TRUE){w<-c(1,0)}	
		
	#Total genetic variance available for selection (Variance of genetic merit/selection goal)

	Vg <-Va+Vm+2*Covam
	Vp <-Va+Vm+Covam+Ve+Vc
	VH <- Vg
	SDH <- sqrt(VH)

	#Covariance of selection goal with phenotype (~index variance)

	covHP <- (Va + (3/2)*Covam + Vm/2)

	#variance of index (bP). Var(Phenotype)=1

	V_I <- covHP								
	SDI <- sqrt(V_I)

	#Accuracy of selection overall

	rHp<-SDI/SDH

	#Overall response to selection with intensity 1

	G_gI <- rHp*SDH

	#regression of genetic values
  
	ba <-(Va + Covam/2)/V_I; bm <- (Vm/2 + Covam)/V_I
	Ga <- ba*SDI; Gm <- bm*SDI
	Resp_a<-0
	Resp_md<-0
	
	DirectResponse[1]<-0
	MaternalResponse[1]<-0
	resultsSUM[1]<-(Ga+Gm)/SDH
	DirectV[1]<-Va
	MaternalV[1]<-Vm
	Covariance[1]<-Covam
	
	#index for results
	counter=2

	#Genetic covariance matrix
	C <- matrix(c(Va, Covam, Covam, Vm), nrow=2)

	#Starting values
	CovaI<-0.00001
	CovmI<-0.00001
	VAs<-0.25*Va
	VAd<-0.25*Va
	Covamd<-Covam/4
	Covams<-Covam/4
	Covam0<-Covam 
	#K's
		
	ks<-calK(ps)
	kd<-calK(pd)
	si<-ks[1]
	di<-kd[1]
	ks<-ks[2]
	kd<-kd[2]
	
	cat("K:s = ", kd, ks , "\n")
	
	Resp_a<-0
	Resp_m<-0
	k=0
	Va0<-Va
	Vm0<-Vm
	VmD<-Vm
	n <- NFS; d <- ND #k=ndam, n=nFS
	raI<-0.00000001; raa<-rmm<-rmI<-raI
	VI<-0.00000001
	
	Pvar=Vp
	

	Covap <- Va + (1/2) * Covam
	CovMdP <- (1/2) * Covam + VmD
	CovaMd <-Covam/2
	K<-(ks+kd)/2
	Va0<-Va*(1-Inbr)
	Vm0<-Vm*(1-Inbr)

	Covam0<-Covam0*(1-Inbr)
	

	
	for(gen in 1:Gen){
		
		
		temp1=NULL
		temp2=NULL
		temp3=NULL
		temp4=NULL
		temp5=NULL
		temp6=NULL
		temp7=NULL
		temp8=NULL
		
				
		D=0
		if(PhenoSelect==FALSE){
		if(IndexSelect==TRUE){
						
			P <- P_matr_index(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
			Pinv <- solve(P)
			if(DirectOnly==TRUE){
				G <- G_matr_index_a(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
				G1<-G[,1]; D=1
				G <- G_matr_index(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
			} else {
				G <- G_matr_index(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
			}
			
		} 
		
		if(DirectOnly==TRUE){
				P <- P_matr_a(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr_a(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
				G1<-G[,1]; D=1
				G <- G_matr(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)				
			} 
		
		if(BLUP==TRUE){
				P <- P_matr(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
			}	
		
		if(independent==TRUE){
				P <- P_matr_independent(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr_independent(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
		
			}
		if(independent_reindeer==TRUE){
				P <- P_matr_independent_reindeer(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr_independent_reindeer(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
		
			}
		
		if(independent2==TRUE){
				P <- P_matr_independent_Wparents(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr_independent_Wparents(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
		
			}
		if(inversedMatingRatio==TRUE){
				P <- P_matr_invMatingR(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr_invMatingR(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
				
		}
		if(Grandsire==TRUE){
				P <- P_matr_Grandsire(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks) 
				Pinv <- solve(P)
				G <- G_matr_Grandsire(d, n, Va0, VAs, VAd, Vm0, VmD, Covam, Covams, Covamd, CovaI, CovmI, kd, ks, VI)
				
		}
		
			if(D==0){
				b<-Pinv %*% G %*% w
				
			} else {
				if (IndexSelect==TRUE){
					b<-c(1, 1/2, 1/4)
					
				}else{
					b<-c(1, 1/2, 1/4,1/2,1/2,0)
				
				}
				
			}
			if(independent==TRUE | inversedMatingRatio==TRUE){
				b[InDbSel==FALSE]<-0
			}
			if(independent2==TRUE){
				b[InDbSel2==FALSE]<-0
			}
			if(independent_reindeer==TRUE){
				b[InDbRein==FALSE]<-0
			}
			
			
			VarI <- t(b) %*% P %*% b
			SD_I <-sqrt(VarI)
			
			if(is.na(SD_I)==TRUE){print(c(n,d, iter));stop()}
			
			rHI <- SD_I/sqrt(Vg)
			#Correlated responses
			
			
			if(DirectOnly==FALSE){
				CovaI <- G[,1]%*%b 
				raI<-CovaI/sqrt(Va*VarI)
				R<-t(G)%*%Pinv%*%G
				
				CovmI <- G[,2]%*%b 
				rmI<-CovmI/sqrt(Vm*VarI)
				raa<-sqrt(R[1,1])/sqrt(Va)
				rmm<-sqrt(R[2,2])/sqrt(Vm)
				ram<-R[1,2]/sqrt(R[1,1]*R[2,2])
			} else {
				CovaI <- G[,1]%*%b 
				covaI <- VarI
				R<-t(G)%*%Pinv%*%G
				raI<-CovaI/sqrt(Va*VarI)
				raa<-sqrt(R[1,1])/sqrt(Va)
				ram<-R[1,2]/sqrt(R[1,1]*R[2,2])
				CovmI<- (G[,2]%*%b)		
				rmI<-CovmI/sqrt(Vm*VarI)
				
				
				rmm<-sqrt(R[2,2])/sqrt(Vm)
				
						
			}
				
			
			
			if(is.na(raa)==TRUE| is.na(rmm)==TRUE){print(c(n,d, iter))}
			covarianses <- c(CovaI,CovmI)
			
			VI <- VarI
			Qresp<-mean(c(si,di))*(covarianses/VI)*SD_I
					
			##add saving variable
			Ga<-sqrt(Va)
			Gm<-sqrt(Vm)
		
		
			
		
			#components affected by seletion (VI before selection)
		
			VmD <- (1-kd*rmI^2)*Vm
			VAs <-0.25*(1-ks*raI^2)*Va; VAd <- 0.25*(1-kd*raI^2)*Va  # between family var comp 
			Covamd <- 0.25*(Covam - (kd*CovaI*CovmI/VI)); Covams <- 0.25*(Covam - (ks*CovaI*CovmI/VI))	
			
			#updating Vm
	
			Vm <- 0.25*(1-ks*rmI^2)*Vm + 0.25*(1-kd*rmI^2)*Vm + Vm0/2
			Covam <- Covams + Covamd + Covam0/2
			Va <- VAs+VAd+Va0/2
			if(DirectOnly==TRUE){
				
			}
			
			Vp <- Va + Vm + Covam + Ve
			
		} else {
			if(gen==1){
				Vp <- as.numeric(P_matr_ph(d,n, Va0,VAs,VAd,Ve, VmD, Covamd,Covam0, Vc,CovaI,CovmI,VI,kd,ks)) 
			}
			
			i <- (si+di)/2
			k <- (ks+kd)/2
			VmD <- Vm
			
			# response in each generatio 0, 1, 2, ...
			# ============================
			covap <- Va + Covam/2
			covmdp <- Covam/2 + Vm

			covmp <- Covam + Vm/2
			# phenotypic selection causing a change in a and m, intensity = 1
			
			Resp_a <- Resp_a + i*covap/sqrt(Vp)
			#cat("Ra",(i*covap/sqrt(Vp)), "\n")
			Resp_m <- Resp_m + i*covmp/sqrt(Vp)

			

			rmdp <- covmdp/sqrt(VmD*Vp)  
			rap <- covap/sqrt(Va*Vp)
			rmp <- covmp/sqrt(Vm*Vp)
			ram <- Covam/sqrt(Va*Vm)


			# var among parents after selection on phenotype at the current gen 0, 1, 2, ...
			Vp <- (1-k)*Vp

			VmD <- (1-kd*rmp^2)*Vm

			# across generations cov(offspring a and dam m)
			
			covamd <- (Covam-kd*rap*rmp*sqrt(Va*Vm))
			Covam <- (Covam -k*rap*rmp*sqrt(Va*Vm))/2 + Covam0/2
			
			# in offspring gen 1, 2, ... (between family + within family var)
			Va <- (1-k*rap^2)*Va/2 + Va0/2
			Vm <- (1-k*rmp^2)*Vm/2 + Vm0/2
			
			#cat("d", Covam)
			# phenotypic var affected by dam's mat var and 2 * cov(offspring a, dam m)
			Vp <- Va + VmD + covamd + Ve


			
			
			Qresp<-c(Resp_a, Resp_m) #responses
		}
		
		
		temp1<-Qresp[1]
		temp2<-raa					#raa<-sqrt(t(G %*% c(1,0))%*%Pinv%*%(G%*%c(1,0)))/sqrt(Va)
		temp3<-Qresp[2]
		temp4<-rmm					#rmm<-sqrt(t(G %*% c(0,1))%*%Pinv%*%(G%*%c(0,1)))/sqrt(Vm)
		temp5<-sum(Qresp)
		temp6<-raI
		temp7<-rmI
		temp8<-VI
		
		DirectV[counter]<-Va
		MaternalV[counter]<-(Vm+VmD)/2
		Covariance[counter]<-Covam
		DirectResponse[counter]<-temp1
		DirectRaa[counter]<-temp2
		MaternalResponse[counter]<-temp3
		MaternalRmm[counter]<-temp4
		resultsSUM[counter]<-temp5
		DirectRaI[counter]<-temp6
		MaternalRmI[counter]<-temp7
		IndexVar[counter]<-temp8
		RAM[counter]<-ram
		if(PhenoSelect==FALSE){
		Pvar[counter]<-Vp 
		} else {
		Pvar[counter]<-Vp
		}
		counter=counter+1
		
	}
	
	Genr<-(0:Gen)
	Results<-cbind(DirectV, MaternalV, Covariance, DirectResponse, DirectRaa, MaternalResponse, MaternalRmm, resultsSUM, DirectRaI, MaternalRmI,IndexVar, Pvar, Genr)
	return(Results)


}
