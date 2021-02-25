


#Discrete generations, 
#rbv
#Object ALL includes pedigree and all sampled values

rm(list=ls())


library(MCMCglmm) 	#rbv, inverseA #csc löytyy
library(MASS)		#csc
library(Matrix) 	#For sparse matrices use Matrix instead of matrix #csc
library(pedigree)	#calcInbreeding 


#3 used in article 1
#4 adds new P ang G

####################################Start of parameters that can be changed#########################################################
#defaults to BLUP
#optional loop for multiple tests, value should be one if only one run is made or FALSE if ignored entirely
multipleTests=10

#Should there be no selection?
PureRandom=FALSE
#Anova for family variance
Anova=FALSE
#BLUP selection in analytical methods (by default stochastic methods revert to BLUP still)
BLUPselection=TRUE
#index selection
IndexSelection=FALSE	#in analytical 
Index_Selection=FALSE#in stochastic
SelectDirectOnly=FALSE #Just direct selection in stochastic (BLUP is the default) Yhden ominaisuuden analyysi.
SDO=FALSE #covariances included in stochastic
WMatr0=FALSE #WW = diag(0)
SelectDirectOnlyAnalytical=FALSE
#Phenotypic selection
PheSelect=FALSE
PheSelectA=FALSE #Remember to adjust weights
#Independent analytical selectors
Independent=FALSE
#with parents
Independent2=FALSE
#for reindeer (only VP and VHS)
Independent_reindeer=FALSE
#inversed mating ratio with analytical selectors. nd is now Males and number of dams is infinite. Only analytical analysis!
InversedMating=FALSE

#Selection on grandsire
Sgrandsire=FALSE
#Save first gen sire raa & rmm
SireAccuracy=FALSE		
#which b are taken into account of the independent selectors. TRUE taken FALSE not
			#P FSmean HSmean
selectors=c(T,T,T)
selectors2=c(T,T,T,T,T,T)
Rselectors=c(T,T)
# variance parameters
Va <- 0.3 #direct effect heritability
Vm <- 0.3 #Indirect effect heritability
rg <- 0.0 #genetic correlation
Covam <- rg * sqrt(Va*Vm)
Vc<-0.0000001

#Family structure
Males=20 								#Number of males
nd=5									#Number of dams per male
nFS=5									#Number of full sibs

#Pedigree and selection
generations=1							#Generations in original pedigree
GenSel=5								#Generations subjected to selection (first generation is reserved for random mating)

#weigths
w<-c(1,1)

##############################################End of parameters#####################################################

###############################################Start of program#####################################################

#Finalizing variance parameters
Ve<-1-(Va+Vm+Covam+Vc)
rg<-Covam/sqrt(Va*Vm)
G<-matrix(c(Va,Covam, Covam, Vm), ncol=2)
E<-matrix(c(Vc,0,0,Ve ), ncol=2)

Alpha<-Ve*(solve(G))

#lambdas
l_a <- Alpha[1,1]; l_m <- Alpha[2,2]; l_am <- Alpha[1,2]

#number of parents
np<-nd*Males+Males


#Ps=(Males)/(0.5*nd*Males*nFS)			#Proportion of selected Males (sires)
#Pd=(nd*Males)/(0.5*nd*Males*nFS)		#Proportion of selected females (dams)

Ps=0.05
Pd=0.75

#For PseudoBLUP, otherwise calK function will fail
if(PureRandom==TRUE){
Ps<-Pd<-.999
}
if(Ps==1){
Ps<-.999
}
if(Pd==1){
Pd<-.999
}

Ft <- 0

#Note: first values in "Analytical[1,]" are of mass selection.
#Objects for saving results
resultsDirect=NULL
resultsMaternal=NULL
resultsCovariance=NULL
resultsInbreeding=NULL
resultsraa=NULL
resultsrmm=NULL
resultsraI=NULL
resultsrmI=NULL
resultsMeanA=NULL
resultsMeanM=NULL
resultsMeasurements=NULL
resultsEstimateCovariance=NULL
familyV=matrix(0, ncol=3, nrow=(generations+GenSel))
#Start of sampling	

if(multipleTests==FALSE){}else{									
for(Mult in 1:multipleTests){

	if(multipleTests>1){
		cat("Test:", Mult, "of", multipleTests, "\n", "...", "\n")
		flush.console()
	}
	

	Pheno=NULL
	temp1=NULL
	temp2=NULL
	raa=NULL
	raI=NULL
	rmm=NULL
	rmI=NULL
	MeanA=NULL
	MeanM=NULL
	
	innerCount=1
			
	females=nd*Males
	np<-Males+females
	ntot<-np+((generations+1)*nFS*females)
		
	#creating original pedigree
	#first 2*np rows are created to give temporary (unrelated) grand parents to parents so that their maternal features can be simulated properly, these rows are removed later. np = number of parents
	
	total<-(females*nFS)
	pedigree<-matrix(NA,ncol=3, nrow=((generations)*(females*nFS)+Males+females+2*np))
	pedigree[,1]<-c(1:nrow(pedigree))
	IDnFS=NULL
	Te=females/Males
		
		
	for(i in 1:(generations)){

		
				
		if(i==1){
			IDmales<-c(1:Males)+2*np
			IDfemales<-c((Males+1):(females+Males))+(2*np)
			pedigree[(2*np+1):(3*np),2]<- 1 : np				#Female ID's are in column 2
			pedigree[(2*np+1):(3*np),3]<-(1+np) : (2*np)
		} else {
			IDmales<-sample(IDnFS, Males)
			IDfemales<-sample(IDnFS[(IDnFS%in%IDmales)==FALSE], females)
		}
				
		IDnFS<-c((3*np+1+((i-1)*total)):(3*np+((i-1)*total)+total))
		matings<-rep(IDmales, Te)
		mates=NULL
					
		mates<-rep(IDmales, each=(Te*nFS))
		NAvecF<-c(rep(IDfemales, each=nFS))
		NAvecM<-c(mates)
		pedigree[c(IDnFS),2]<-NAvecF
		pedigree[c(IDnFS),3]<-NAvecM
				
	}
		
		
		
	ID<-pedigree[,1] 
	Pheno<-rbv(pedigree, G)
		
	x<-unique(pedigree[,2])[-1]
	temp=matrix(rep(0, 3*length(ID)), ncol=3)
	
	#Simulating environmental effects and setting the visible maternal effects
	N=1
	for(xx in 1:length(x)){
	
		if(x[xx] > np){N<-nFS}
		temp[pedigree[,2]%in%x[xx],1]<-rnorm(N, 0, sqrt(Ve))
		temp[pedigree[,2]%in%x[xx],2]<-rep(rnorm(1, 0,  sqrt(Vc)),N) ##family Vc
		temp[pedigree[,2]%in%x[xx],3]<-Pheno[x[xx],2] #Maternal effects
	}
		
	Pheno<-cbind(Pheno,temp)[-1:-(2*np),]
	pedigree<-pedigree[-1:-(2*np),]
	ALL2<-pedigree
	pedigree<-pedigree-(2*np)
	
	pedigree[1:np,2:3]<-NA
	Pedigree<-pedigree
	measurements<-as.vector(rowSums(Pheno[,c(1,3,4,5)]))
	
	inbreeding<-calcInbreeding(as.data.frame(pedigree))
	Inbreeding<-mean(tail(inbreeding, total))
	ALL<-cbind(pedigree, Pheno, measurements,inbreeding)
	ALL2<-cbind(ALL2,Pheno, measurements,inbreeding)
	colnames(ALL)<-c("ID", "DAM", "SIRE", "Ga", "Gm", "Pe", "Pc", "Pm", "measurements", "inbreeding")
	rownames(ALL)<-c()
	colnames(ALL2)<-c("ID", "DAM", "SIRE", "Ga", "Gm", "Pe", "Pc", "Pm", "measurements", "inbreeding")
	rownames(ALL2)<-c()
	MeanA<-mean(ALL[,4])
	MeanM<-mean(ALL[,5])
	#variance components
	VA<-Va; VM<-Vm; COV<-Covam; VP<-var(measurements)
	TempOffs=NULL
		
	#Save empirical variances
	VA<-c(VA, var(tail(ALL[,4], total)))
	VM<-c(VM, var(tail(ALL[,5], total)))
	COV<-c(COV, cov(tail(ALL[,4], total),tail(ALL[,5], total)))
		
	for(SelGen in 1:GenSel){
			
		#cat("Va:", var(tail(ALL[,4], total)), " Vm:", var(tail(ALL[,5], total)), 
		#	" Covam:", cov(tail(ALL[,4], total), tail(ALL[,5], total)),"\n")
		cat("Generation under selection: ", SelGen, "\n")
		#cat("Variance components after a round of selection:", "\n")
		flush.console()
			
		ntot<-np+((generations+SelGen-1)*nFS*females)
		if(Index_Selection==TRUE){
						
			Dim<-dim(pedigree)[1]
			pedigree[ -(Dim-total+1):-Dim, 2:3]<-NA
		
		}
			#ordering the selection candidates
			# 1. to Males and females with p(male)=0.5 and p(female)=1-p(male)
			calfM<-sample(1:total, total/2)
			calfF<-c(1:total)[-calfM]
			if(SelGen==1){
				ssss<-calfM
				ffff<-calfF
			}
		if(PheSelect==FALSE){
			
			#incidence matrices
			#observation vector for direct effects 
			y <- c(rep(0,np),rep(1,(ntot-np)))
			
			Z<-Matrix(0, nrow=ntot,ncol=ntot, sparse=TRUE)
			diag(Z) <- y
			
			
			Z<-Matrix(Z, sparse=TRUE)
						
			ZZ<-Z 		#=Z
						
			#maternal effects
			Xcoord<-pedigree[(np+1):ntot,2]
			#Ycoord<-pedigree[(np+1):ntot,1]
			W <- Matrix(0,nrow=ntot,ncol=ntot ,sparse=TRUE)
			
			if(WMatr0==FALSE){
			if(Index_Selection==TRUE){  #redundant
				for(i in (Dim-total+1):Dim){
					W[pedigree[i,1],pedigree[i,2]]<-1
				}
			} else {
				for(i in c((np+1):ntot)){
					W[pedigree[i,1],pedigree[i,2]]<-1
				}
			}
			}
			
			
			WW <- crossprod(W,W)
		
			# direct / maternal effects
			ZW = crossprod(Z, W)
			
			WZ = crossprod(W,Z)
			
			
			invA <- inverseA(pedigree)
			Ainv <- invA$Ainv
			
			
			ObsErA<-(t(Z)) %*% measurements
			ObsErM<-(t(W)) %*% measurements
			
			if(SelectDirectOnly==FALSE){
				ObsEr<-rbind(ObsErA,ObsErM)
				
				flush.console()
				# mixed model equations
								
				MME11 <- ZZ + Ainv*l_a
				MME12 <- ZW + Ainv*l_am
				MME21 <- WZ + Ainv*l_am
				MME22 <- WW + Ainv*l_m
				
				Dimension <- dim(MME11)[1]
				MME <- Matrix(0, nrow = (Dimension*2), ncol= (Dimension*2), sparse=TRUE)
				
							
				MME[1:Dimension, 1:Dimension] <- MME11					; 	MME[1:Dimension, (Dimension+1):(2*Dimension)] <- MME12
				MME[(Dimension+1):(2*Dimension), 1:Dimension] <- MME21	; 	MME[(Dimension+1):(2*Dimension), (Dimension+1):(2*Dimension)] <- MME22
				
				
				#memory purge
				ZZ<-ZW<-WZ<-WW<-NULL
				
				#calculating predictions
				
				predictions <-solve(MME, ObsEr)
				
				#Taking out the separate predictions for direct and maternal effects
				a<-length(predictions)
				predA<-predictions[1:(a/2)]
				predM<-predictions[((a/2)+1):(a)]
				COR<-cor(predA,predM)
				if(SDO==TRUE){predM<-rep(0, length(predA))}
			} else {
				
				
				MME<-(Ainv)
				MME <- Matrix(MME, sparse=TRUE)
				
				#calculating predictions
				predictions <- solve(MME, ObsErA)
				predA<-as.vector(predictions)			
				predM<-rep(0, length(predA))
				COR=0
			}
			
			
			pI<-predA*w[1]+predM*w[2]
			
			#Forming an index
		
		
			#calculating correlations
			
			raa<-c(raa,cor(tail(predA,total), tail(ALL[,4],total)))
			raI<-c(raI,cor(tail(pI,total), tail(ALL[,4],total)))
		
			
			rmm<-c(rmm,cor(tail(predM,total), tail(ALL[,5],total)))
			rmI<-c(rmI,cor(tail(pI,total), tail(ALL[,5],total)))
			#Calculating empirical index variance and empirical covariances between effects and index
			VpI<-var(tail(pI,total))
			caI<-cov(tail(predA,total),tail(Pheno[,1],total))
			cmI<-cov(tail(predictions[-1:-np],total),tail(Pheno[,2],total))
		
			
			
		
			if( SelGen!= 1  & PureRandom == FALSE){
				# 2. Order from best to worst by sex and select the breeding population
				bestM<-order(tail(pI, total)[calfM], decreasing=TRUE) #have to distinct the sexes and pick accordingly
				bestM<-bestM[1:Males]
				bestF<-order(tail(pI, total)[calfF], decreasing=TRUE)
				bestF<-bestF[1:(nd*Males)]
			
				sires<-tail(ALL, total)[calfM,][bestM,] #These are now ordered from best to worst
				dams<-tail(ALL, total)[calfF,][bestF,]
			
			
				sires<-sires[sample(1:nrow(sires),),] #Now they are randomized again. (Without this step there seemed to be an additional slight positive trend for mean over generations)
				dams<-dams[sample(1:nrow(dams),),]

			
				#random mating between selected population
				sT<-1:Males
				dT<-1:females
			} else {
				sires<-tail(ALL, total)[calfM,]
				dams<-tail(ALL, total)[calfF,]
				#random mating between randomly selected animals
				sT<-sample(1:nrow(sires), Males, replace=FALSE)
				dT<-sample(1:nrow(dams), females, replace=FALSE)
		
			}	
		} else {
			#PheSelect
			#calculating correlations, these are now just sqrt(heritabilities) but called raa, raI etc for the stability of the code
			#raa<-c(raa,cor(predA, ALL[,4]))
			raa<-c(raa,cor(tail(measurements,total), tail(ALL[,4],total)))
			raI<-c(raI,cor(tail(measurements,total), tail(ALL[,4],total)))
		
			#rmm<-c(rmm,cor(predM, ALL[,5]))
			rmm<-c(rmm,cor(tail(measurements,total), tail(ALL[,5],total)))
			rmI<-c(rmI,cor(tail(measurements,total), tail(ALL[,5],total)))
			bestM<-order(tail(measurements, total)[calfM], decreasing=TRUE)
			bestM<-bestM[1:Males]
			bestF<-order(tail(measurements, total)[calfF], decreasing=TRUE)
			bestF<-bestF[1:(nd*Males)]
			
			sires<-tail(ALL, total)[calfM,][bestM,] #These are now ordered from best to worst
			dams<-tail(ALL, total)[calfF,][bestF,]
			###################################################################MeanM<-c(MeanM, mean(c(sires[,4], dams[,4]))) #########test remove
			sires<-sires[sample(1:nrow(sires),),] #Now they are randomized again. (Without this step there seemed to be an additional slight positive trend for mean over generations)
			dams<-dams[sample(1:nrow(dams),),]
			sT<-1:Males
			dT<-1:females
		}
		
		
		
		for(k in 1:Males){ 
				
			Asegr=NULL
			Msegr=NULL
			ratio<-(females/Males)
			#Pick male and females
			if(length(sT)>1 & length(dT)>1){
			snd<-sample(sT,1, replace=FALSE)
			fnd<-sample(dT, ratio, replace=FALSE)
			} else {
			snd<-sT
			fnd<-dT
			}
			sT<-sT[sT!=snd]
			dT<-dT[!(dT %in% fnd)]
				
			#Calculate values for effects
			IDs<-(nrow(ALL)+1):(nrow(ALL)+(nFS*ratio))					#animal ID
			PindexS<-rep(sires[snd,1],(nFS*ratio))						#sire ID
			PindexF<-as.vector(rep(dams[fnd,1],each=nFS))				#Dam IDs
			
			pedigree<-rbind(pedigree, cbind(IDs, PindexF, PindexS))		#Extending pedigree
			Pedigree<-rbind(Pedigree, cbind(IDs, PindexF, PindexS))
			
			dat<-as.data.frame(pedigree)
			inbreeding<-calcInbreeding(dat)
			Fs<-tail(inbreeding, length(IDs))
				
			Aherited<-as.vector(0.5*sires[snd,4]+0.5*dams[fnd,4])		#0.5As + 0.5Ad
			Aherited<-rep(Aherited, each=nFS)
			Mherited<-as.vector(0.5*sires[snd,5]+0.5*dams[fnd,5])		#0.5Ms +0.5Md
			Mherited<-rep(Mherited, each=nFS)
			Mobs<-rep(ALL[unique(dams[fnd,1]),5], each=nFS)				#m observed	
				
			VaSegr<-0.5*(1-Fs)*Va
			VmSegr<-0.5*(1-Fs)*Vm
			CovSegr<-0.5*(1-Fs)*Covam
			
			
				
			x=NULL
				
			for(ze in 1:length(IDs)){
				a<-VaSegr[ze]
				m<-VmSegr[ze]
				Sc<-CovSegr[ze]
				Gi<-matrix(c(a,Sc,Sc,m),ncol=2)
				x<-rbind(x, mvrnorm(n=1, mu=c(0, 0), Sigma=Gi))
			}
							
			Asegr<-x[,1]												#Asegr (0.5 A0)
			Msegr<-x[,2]
			
			Common<-rep(rnorm(ratio, mean=0, sd=sqrt(Vc)), each=nFS)	#Common environment + (maternal segregation)
			Env<-rnorm(ratio*nFS, mean=0, sd=sqrt(Ve))					#E
			
			Obs<-Aherited+Asegr+Mobs+Common+Env				#Sums
				
				
			#Extend pedigree and effects table as well as the measurements vector
			ALL<-rbind(ALL, cbind(IDs, PindexF, PindexS, (Aherited+Asegr), (Mherited+Msegr), Env, Common, Mobs, Obs,Fs)) 
				
			measurements<-c(measurements, Obs)
			
		}
			dat<-as.data.frame(Pedigree)
			
			ALL[,10]<-calcInbreeding(dat)
			
		#Saving variance objects
		VA<-c(VA, var(tail(ALL[,4], total)))
		VM<-c(VM, var(tail(ALL[,5], total)))
		VP<-c(VP, var(tail(measurements, total)))
		COV<-c(COV, cov(tail(ALL[,4], total),tail(ALL[,5], total)))
		MeanA<-c(MeanA, mean(tail(ALL[,4], total)))
		MeanM<-c(MeanM, mean(tail(ALL[,5], total)))
		Inbreeding<-c(Inbreeding, mean(tail(ALL[,10], total)))
		if(SelGen==GenSel){
			#cat("Va:", var(tail(ALL[,4], total)), " Vm:", var(tail(ALL[,5], total)), 
			#" Covam:", cov(tail(ALL[,4], total), tail(ALL[,5], total)),"\n")
		}
		
	}
	
	if(nd>1 & Anova==TRUE){
		
			Temp<-ALL[-1:-np,]
			Temp[,2:3]<-Temp[,2:3]+(2*np)
			ALL2<-rbind(ALL2,Temp)
			if(Mult==1){
			familyV<-Vfamilies(dat=ALL, depth=(generations+GenSel), NP=np, N=total)
			familyV<-familyV/multipleTests
			} else {
			Temp<-Vfamilies(dat=ALL, depth=(generations+GenSel), NP=np, N=total)
			Temp<-Temp/multipleTests
			familyV <- familyV + Temp
			}
		}
		

	#VA	
	#VM
	
	
	resultsDirect<-rbind(resultsDirect, VA)
	resultsMaternal<-rbind(resultsMaternal,VM)
	resultsCovariance<-rbind(resultsCovariance,COV)
	resultsInbreeding<-rbind(resultsInbreeding, Inbreeding)
	resultsraa<-rbind(resultsraa,raa)
	resultsraI<-rbind(resultsraI,raI)
	resultsrmm<-rbind(resultsrmm,rmm)
	resultsrmI<-rbind(resultsrmI,rmI)
	resultsMeanA<-rbind(resultsMeanA, MeanA)
	resultsMeanM<-rbind(resultsMeanM, MeanM)
	resultsMeasurements<-rbind(resultsMeasurements, VP)
	
}	

Inbreeding<-colMeans(resultsInbreeding)
Ft<-Inbreeding[length(Inbreeding)-1]
}




cat("Inbreeding overall = ", Ft, "\n")
Results<-cbind(resultsMeanA, resultsMeanM)
write.table(Results, file="filename.txt", quote=FALSE, sep="\t", eol="\n")


#Additional code for ease of use (first place in responce vector is reserved for empirical starting value for the response, secon is reserved for first generation of random mating)

colMeans(resultsMeanA) #direct effect response
colMeans(resultsMeanM) #indirect effect response
colMeans(resultsraa)	#correlation between estimated and true direct effect
colMeans(resultsrmm)	#correlation between estimeted and true indirect effect	


















