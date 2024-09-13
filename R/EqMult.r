



# function to minimize for MIRF or MTRF method
obj<-function(par,T,ni,ab,aj1T,bj1T,cj1T,met,itmp,wt,D=D,base){
	A<-rep(1,T)
	A[-base]<-par[1:(T-1)]
	B<-rep(0,T)
	B[-base]<-par[T:(2*T-2)]
	if (itmp==1) A[1:T]<-1
	bj1Ts<-bj1T
	for (t in 1:T) bj1Ts[,t]<-(bj1Ts[,t]*A[t]+B[t])
	aj1Ts<-aj1T
	for (t in 1:T) aj1Ts[,t]<-aj1Ts[,t]/A[t]
	bj<-rowMeans(bj1Ts,na.rm=T)
	aj<-rowMeans(aj1Ts,na.rm=T)

	f<-0
	for (t in 1:T) {
		bjts<-(bj-B[t])/A[t]
		ajts<-aj*A[t]
		ajt<-aj1T[,t]
		bjt<-bj1T[,t]
		cjt<-cj1T[,t]
		Pt<-matrix(NA,length(ab),ni)
		for (i in 1:ni)
			Pt[,i]<-equateIRT:::irtp1(ab,diff=bjt[i],discr=ajt[i],guess=cjt[i],D=D)
		Pts<-matrix(NA,length(ab),ni)
		for (i in 1:ni)
			Pts[,i]<-equateIRT:::irtp1(ab,diff=bjts[i],discr=ajts[i],guess=cjt[i],D=D)
		Pts[is.na(Pt)]<-NA
		if (met=="irf") f<-f+0.5*sum((Pt-Pts)^2*wt,na.rm=T)
		if (met=="trf") f<-f+0.5*sum(((rowSums(Pt,na.rm=T)-rowSums(Pts,na.rm=T))^2)*wt)
	}

	return(f)
}

# numerical first derivative of obj with respect to the item parameters
# used for finding the numerical second derivative of obj
# with respect to the equating coefficients and the item parameters
Sgamma<-function(x,par,T,ni,ab,tabl,nummet,itmp,wt,D=D,base){
	tabl$gamma<-x[rownames(tabl)]
	tab<-reshape(tabl,direction="wide",v.names="gamma",timevar="time", idvar = "itms")
	bj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Dffclt",][,-1])
	if (itmp>1) aj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Dscrmn",][,-1])
	else aj1T<-matrix(1,ni,T)
	if (itmp==3) cj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Gussng",][,-1])
	else cj1T<-matrix(0,ni,T)
	grad(func=objectivefzRcpp,x=par,T=T,ab=ab,wt=wt,aj1T=aj1T,bj1T=bj1T,cj1T=cj1T,nummet=nummet,itmp=itmp,D=D,base=base) 
}



# this functions returns the matrix of second derivatives of obj
# (at the minimum) with respect to the equating coefficients and 
# the item parameters
# calls Rcpp
pABgammaR2C<-function(par,T,ab,wt,aj1T,bj1T,cj1T,nummet,itmp,D,base)
{
	nb<-sum(!is.na(bj1T))
	posnomi<-bj1T
	posnomi[!is.na(posnomi)]<-1:nb
	tabnomia<-outer(rownames(aj1T),1:T,paste,sep=".")
	tabnomia[is.na(aj1T)]<-NA
	tabnomib<-outer(rownames(bj1T),1:T,paste,sep=".")
	tabnomib[is.na(bj1T)]<-NA
	if (itmp==3) { tabnomic<-outer(rownames(cj1T),1:T,paste,sep=".")
		tabnomic[is.na(cj1T)]<-NA}
	nomia<-tabnomia[!is.na(tabnomia)]
	nomib<-tabnomib[!is.na(tabnomib)]
	if (itmp==3) nomic<-tabnomic[!is.na(tabnomic)]

	tt<-partialABgammaRcpp(par,T,ab,wt,aj1T,bj1T,cj1T,nummet,itmp,D,base,nb,posnomi)

	out_a<-tt[[1]]
	out_b<-tt[[2]]
	if(itmp==3) out_c<-tt[[3]]
	else out_c<-NULL
	colnames(out_a)<-nomia
	colnames(out_b)<-nomib
	if (itmp==3) colnames(out_c)<-nomic

	out<-cbind(out_c,out_a,out_b)
	out<-out[-c(base,base+T),]
	return(out)
}



# this function is used to compute the numerical derivatives of a_j^* and b_j^* 
# with respect to the item parameter estimates
# as input requires the item parameter estimates (gamma)
# as output returns a_j^* and b_j^* 
gamma2ab<-function(x,par,T,ni,ab,tabl,nummet,itmp,wt,D=D,base,ini){
	#print(x)
	#Ts<-paste("T",1:T,sep="")
	tabl$gamma<-x[rownames(tabl)]
	tab<-reshape(tabl,direction="wide",v.names="gamma",timevar="time", idvar = "itms")
	bj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Dffclt",][,-1])
	if (itmp>1) aj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Dscrmn",][,-1])
	else aj1T<-matrix(1,ni,T)
	if (itmp==3) cj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Gussng",][,-1])
	else cj1T<-matrix(0,ni,T)
	out<-nlminb(start=ini,objective=objectivefzRcpp,T=T,ab=ab,wt=wt,aj1T=aj1T,bj1T=bj1T,cj1T=cj1T,nummet=nummet,itmp=itmp,D=1,base=base,control=list(iter.max=100000,eval.max=10000,trace=FALSE))
	par<-out$par
	As<-rep(1,T)
	#names(As)<-Ts
	As[-base]<-par[1:(T-1)]
	Bs<-rep(0,T)
	#names(Bs)<-Ts
	Bs[-base]<-par[T:(2*T-2)]
	bj1Ts<-bj1T
	for (t in 1:T) bj1Ts[,t]<-bj1Ts[,t]*As[t]+Bs[t]
	aj1Ts<-aj1T
	for (t in 1:T) aj1Ts[,t]<-aj1Ts[,t]/As[t]
	bj<-rowMeans(bj1Ts,na.rm=T)
	aj<-rowMeans(aj1Ts,na.rm=T)
	#names(bj)<-tab$itms[substr(tab$itms,1,6)=="Dffclt"]
	#names(aj)<-tab$itms[substr(tab$itms,1,6)=="Dscrmn"]
	return(c(aj,bj))
}



multiec_irf<-function(mods,base,nq=30,method,se,start,iter.max,trace=FALSE){
	cat("Computation of equating coefficients  .  .  .  . \n")
	if (inherits(start,"mlteqc")) start<-c(start$A[-1],start$B[-1])
	itms<-c()
	for (k in 1:length(mods)) itms<-c(itms,names(mods[[k]]$coef))
	itms<-sort(unique(itms))
	T<-length(mods)
	#Ts<-paste("T",1:T,sep="")
	modsnames<-names(mods)
	tab<-data.frame(itms=itms)
	for (k in 1:length(mods)) 
		tab<-merge(tab,mods[[k]]$coef,by.x="itms",by.y=0,all=T,suffixes=c(k-1,k))
	colnames(tab)[-1]<-modsnames
	tab$itms<-as.character(tab$itms)
	rownames(tab)<-tab$itms
	sel_administered2plus<-rowSums(!is.na(tab[,-1]))>=2
	items_administered2plus<-tab$itms[sel_administered2plus]
	tab_restr<-tab[sel_administered2plus,]
	
	itmp<-2
	if (sum(substr(tab$itms,1,6)=="Dscrmn")==0) itmp=1
	if (sum(substr(tab$itms,1,6)=="Gussng")>0) itmp=3

	if (is.null(start)) ini<-c(rep(1,T-1),rep(0,T-1))
	if (!is.null(start)) ini<-start
	names(ini)<-c(modsnames[-base],modsnames[-base])
	
	gq<-gauss.quad.prob(nq,dist="normal")
	ab<-gq$nodes
	wt<-gq$weights

	ni<-nrow(tab)/itmp
	ni_restr<-nrow(tab_restr)/itmp
	bj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Dffclt",][,-1])
	if (itmp>1) aj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Dscrmn",][,-1])
	if (itmp==1) {aj1T<-matrix(1,ni,T); aj1T[is.na(bj1T)]<-NA}
	if (itmp==3) cj1T<-as.matrix(tab[substr(tab$itms,1,6)=="Gussng",][,-1])
	if (itmp<3) cj1T<-matrix(0,ni,T)
	if (method=="irf") nummet<-1
	if (method=="trf") nummet<-2
	bj1T_restr<-bj1T[rownames(bj1T)%in%items_administered2plus,]
	if (itmp>1) aj1T_restr<-aj1T[rownames(aj1T)%in%items_administered2plus,]
	if (itmp==1) {aj1T_restr<-matrix(1,ni_restr,T); aj1T[is.na(bj1T)]<-NA}
	if (itmp==3) cj1T_restr<-cj1T[rownames(cj1T)%in%items_administered2plus,]
	if (itmp<3) cj1T_restr<-matrix(0,ni_restr,T)
	# the following uses Rcpp function
	out<-nlminb(start=ini,objective=objectivefzRcpp,gradient=gradRcpp,hessian=hessRcpp,T=T,ab=ab,wt=wt,aj1T=aj1T_restr,bj1T=bj1T_restr,cj1T=cj1T_restr,nummet=nummet,itmp=itmp,D=1,base=base,control=list(iter.max=iter.max,trace=trace))
	par<-out$par
	As<-rep(1,T)
	names(As)<-modsnames
	if (itmp>1) As[-base]<-par[1:(T-1)]
	Bs<-rep(0,T)
	names(Bs)<-modsnames
	Bs[-base]<-par[T:(2*T-2)]
	
	bj1Ts<-bj1T
	for (t in 1:T) bj1Ts[,t]<-bj1Ts[,t]*As[t]+Bs[t]
	aj1Ts<-aj1T
	for (t in 1:T) aj1Ts[,t]<-aj1Ts[,t]/As[t]
	bj<-rowMeans(bj1Ts,na.rm=T)
	aj<-rowMeans(aj1Ts,na.rm=T)
	names(bj)<-tab$itms[substr(tab$itms,1,6)=="Dffclt"]
	names(aj)<-tab$itms[substr(tab$itms,1,6)=="Dscrmn"]
	rownames(aj1T)<-tab$itms[substr(tab$itms,1,6)=="Dscrmn"]
	rownames(bj1T)<-tab$itms[substr(tab$itms,1,6)=="Dffclt"]
	rownames(cj1T)<-tab$itms[substr(tab$itms,1,6)=="Gussng"]
	
	if (se) {
		cat("Computation of standard errors ")
		# h1 is the matrix of second derivatives at the minimum (hessian)
		# the following uses Rcpp function
		h1<-hessRcpp(par=out$par,T=T,ab=ab,wt=wt,aj1T=aj1T_restr,bj1T=bj1T_restr,cj1T=cj1T_restr,nummet=nummet,itmp=itmp,D=1,base=base)
		colnames(h1)<-c(rep("A",T-1),rep("B",T-1))
		cat(" . ")

		# pfABg is the matrix of second derivatives with respect to the equating coefficients and the item parameters
		# the following uses Rcpp function
		pfABg<-pABgammaR2C(par=out$par,T=T,ab=ab,wt=wt,aj1T=aj1T_restr,bj1T=bj1T_restr,cj1T=cj1T_restr,nummet=nummet,itmp=itmp,D=1,base=base)
		if (itmp==1) {pfABg[1:(T-1),]<-0; pfABg<-pfABg[,grep("Dffclt",colnames(pfABg))]}
		cat(" . ")

		if (itmp==1) {h1<-h1[-(1:(T-1)),-(1:(T-1))]; pfABg<-pfABg[-(1:(T-1)),]}
		dAB_gamma<--solve(h1)%*%pfABg
		dAB_gamma<-dAB_gamma[,colSums(dAB_gamma)!=0]
		cat(" . ")

		sel<-colnames(dAB_gamma)
		vars<-lapply(mods,FUN=function(x) x$var)
		for (i in 1:length(vars)) rownames(vars[[i]])<-colnames(vars[[i]])<-paste(rownames(vars[[i]]),i,sep=".")
		VarAll<-VarExtRcpp(vars)
		VarAllNames<-c()
		for (i in 1:length(vars)) VarAllNames<-c(VarAllNames,rownames(vars[[i]]))
		rownames(VarAll)<-colnames(VarAll)<-VarAllNames
		#VarAll<-VarExt(mods)
		varAB<-dAB_gamma%*%VarAll[sel,sel]%*%t(dAB_gamma)
		seAB<-diag(varAB)^0.5
		seA<-rep(0,T)
		names(seA)<-modsnames
		if (itmp>1) seA[-base]<-seAB[1:(T-1)]
		seB<-rep(0,T)
		names(seB)<-modsnames
		if (itmp>1) seB[-base]<-seAB[T:(2*T-2)]
		if (itmp==1) seB[-base]<-seAB[1:(T-1)]
		
		# standard errors of synthetic item parameters
		
		invT<-as.matrix(rowSums(!is.na(bj1T))) # -> u_j
		invT<-invT[,c(rep(1,T))]
		colnames(invT)<-modsnames

		if (itmp>1) {
			pajA<-as.matrix(aj1T)
			for (t in 1:T) pajA[,t]<-aj1T[,t]/As[t]^2
			pajA[is.na(pajA)]<-0
			pajA<- -pajA/invT
			pajA<-as.matrix(pajA[,-base])
			
			tablong<-reshape(tab,direction="long",varying = list(2:6),idvar="itms",v.names="value")
			tablong<-tablong[!is.na(tablong$value),]
			itms_t<-rownames(tablong)
			dAB_gamma_all<-matrix(0,(T-1)*2,length(itms_t))
			colnames(dAB_gamma_all)<-itms_t
			rownames(dAB_gamma_all)<-rownames(dAB_gamma)
			dAB_gamma_all[,colnames(dAB_gamma)]<-dAB_gamma
			
			grAa<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="A"),select=grepl("^Dscrmn",colnames(dAB_gamma_all))) # derivatives of equating coefficients A with respect to the discrimination item parameters
			pajgamma_a<-pajA%*%grAa # this is the result of the sum in Equations (91) and (92) multiplied by -1/u_j
			
			colnames_spl<-strsplit(colnames(grAa),split=".",fixed=T)
			colnames_no_t<-sapply(colnames_spl,FUN=function(x) paste(x[1],x[2],sep="."))
			sel_sameitem<-outer(rownames(aj1T),colnames_no_t,FUN="==")
			
			whicht<-as.numeric(as.character(sapply(colnames_spl,FUN=function(x) x[3])))
			
			pajgamma_a[sel_sameitem]<-pajgamma_a[sel_sameitem]+(matrix(1/invT[,1])%*%(1/As[whicht]))[sel_sameitem] # these are the derivatives in Equations (91) and (92)

			grAb<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="A"),select=grepl("^Dffclt",colnames(dAB_gamma_all))) # derivatives of equating coefficients A with respect to the difficulty item parameters
			pajgamma_b<-pajA%*%grAb # this is the result of the sum in Equation (93) multiplied by -1/u_j

			grAc<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="A"),select=grepl("^Gussng",colnames(dAB_gamma_all))) # derivatives of equating coefficients A with respect to the difficulty item parameters
			if (itmp==3) pajgamma_c<-pajA%*%grAc # similar to Equation (92), derivative of a_j^* with respect to c_it
			else pajgamma_c<-c()

			pajgamma<-cbind(pajgamma_a,pajgamma_b,pajgamma_c)
			
			sel<-colnames(pajgamma)
			seaj<-diag(pajgamma%*%VarAll[sel,sel]%*%t(pajgamma))^0.5
			seaj<-seaj[names(aj)]
		}
		if (itmp==1) aj<-seaj<-NULL
		
		cat(" . \n")
		
		pbjA<-as.matrix(bj1T)
		pbjA[is.na(pbjA)]<-0
		pbjA<- pbjA/invT
		pbjA<-as.matrix(pbjA[,-base])
		
		pbjB<-pbjA
		pbjB[pbjB!=0]<-1/invT[,-base][pbjB!=0]
		
		grAa<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="A"),select=grepl("^Dscrmn",colnames(dAB_gamma_all))) # derivatives of equating coefficients A with respect to the discrimination item parameters
		grBa<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="B"),select=grepl("^Dscrmn",colnames(dAB_gamma_all))) # derivatives of equating coefficients B with respect to the discrimination item parameters
		if (itmp>1) pbjgamma_a<-pbjA%*%grAa+pbjB%*%grBa # Equation (94)
		if (itmp==1) pbjgamma_a<-pbjB%*%grBa # Equation (94)
		
		grAb<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="A"),select=grepl("^Dffclt",colnames(dAB_gamma_all))) # derivatives of equating coefficients A with respect to the difficulty item parameters
		grBb<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="B"),select=grepl("^Dffclt",colnames(dAB_gamma_all))) # derivatives of equating coefficients B with respect to the difficulty item parameters
		if (itmp>1) pbjgamma_b<-pbjA%*%grAb+pbjB%*%grBb # this is the result of the sum in Equations (95) and (96) multiplied by -1/u_j
		if (itmp==1) pbjgamma_b<-pbjB%*%grBb # this is the result of the sum in Equations (95) and (96) multiplied by -1/u_j

		grAc<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="A"),select=grepl("^Gussng",colnames(dAB_gamma_all))) # derivatives of equating coefficients A with respect to the guessing item parameters
		grBc<-subset(dAB_gamma_all,subset=(rownames(dAB_gamma_all)=="B"),select=grepl("^Gussng",colnames(dAB_gamma_all))) # derivatives of equating coefficients B with respect to the guessing item parameters
		if (itmp==3) pbjgamma_c<-pbjA%*%grAc+pbjB%*%grBc # similar to Equation (94), derivative of b_j^* with respect to c_it
		else pbjgamma_c<-c()

		colnames_spl<-strsplit(colnames(grBb),split=".",fixed=T)
		colnames_no_t<-sapply(colnames_spl,FUN=function(x) paste(x[1],x[2],sep="."))
		sel_sameitem<-outer(rownames(bj1T_restr),colnames_no_t,FUN="==")
		
		whicht<-as.numeric(as.character(sapply(colnames_spl,FUN=function(x) x[3])))
		
		pbjgamma_b[sel_sameitem]<-pbjgamma_b[sel_sameitem]+(matrix(1/invT[,1])%*%As[whicht])[sel_sameitem]
		
		pbjgamma<-cbind(pbjgamma_a,pbjgamma_b,pbjgamma_c)
		
		sel<-colnames(pbjgamma)
		sebj<-diag(pbjgamma%*%VarAll[sel,sel]%*%t(pbjgamma))^0.5
		sebj<-sebj[names(bj)]
		
		partial<-t(dAB_gamma) # derivatives of A and B equating coefficients with respect to the item parameters
		if (itmp==1) partial<-cbind(matrix(0,nrow(partial),T-1),partial) 
		namesAB<-c(paste("A",(1:T)[-base],sep="."),paste("B",(1:T)[-base],sep="."))
		colnames(partial)<-namesAB
		if (itmp==1)
		{
		  varAB1<-matrix(0,2*T-2,2*T-2)
		  varAB1[T:(2*T-2),T:(2*T-2)]<-varAB
		  varAB<-varAB1
		}
		rownames(varAB)<-colnames(varAB)<-namesAB
	}
	else {
		seA<-rep(NA,T)
		seB<-rep(NA,T)
		varAB<-matrix(NA,T,T)
		seaj<-NULL
		sebj<-NULL
		vars<-NULL
		partial<-NULL
	}
	if (itmp==1) tabs<-bj1Ts
	if (itmp==2) tabs<-rbind(bj1Ts,aj1Ts)
	if (itmp==3) tabs<-rbind(bj1Ts,aj1Ts,cj1T)
	colnames(tabs)<-paste(colnames(tabs),modsnames[base],sep=".as.")
	tabs<-tabs[,-base]
	tab<-cbind(tab,tabs)
	colnames(tab)[1]<-"Item"
	out<-list(A=As,B=Bs,se.A=seA,se.B=seB,varAB=varAB,as=aj,bs=bj,se.as=seaj,se.bs=sebj,tab=tab,varFull=vars,partial=partial,itmp=itmp,method=method,basename=modsnames[base],convergence=out$convergence)
	class(out)<-"mlteqc"
	return(out)
}






ipf4der<-function(gamma,itms,t,base,aj1T) {
	names(gamma)<-itms
	for (tt in 1:max(t)) aj1T[,paste("gamma",tt,sep=".")]<-gamma[t==tt][rownames(aj1T)]
	mmout<-ipfRcpp(as.matrix(aj1T[,-1]),base,0.00001)
	As<-mmout[[1]]
	return(As)
}



# this function is used to calculate numerical derivatives
# of a*, A, b*, B with respect to the discrimination and difficulty parameters
der_asAbsB_ab<-function(gamma,tab,base,method1,T,ni)
{
	tab$gamma<-gamma
	itms<-tab$itempar
	#Ts<-paste("T",1:T,sep="") # label of administrations: T1, T2, ...

	tabDscrmn<-tab[substr(tab$itempar,1,6)=="Dscrmn",] # table with discrimination parameters
	tabDffclt<-tab[substr(tab$itempar,1,6)=="Dffclt",] # table with difficulty parameters
	
	# MULTIPLE MEAN-GEOMETRIC MEAN (performed anyway)
	# step 1: estimation of A equating coefficients
	if (nrow(tabDscrmn)>0) {
		if (any(tabDscrmn$gamma<0)) warning("Negative discrimination parameter. Converted to positive value.")
		tabDscrmn$gamma<-abs(tabDscrmn$gamma) # negative discriminations converted to positive values
		tabDscrmn$t<-as.factor(tabDscrmn$t)
		tabDscrmn$t <- relevel(tabDscrmn$t, ref = base) #set the A equating coefficient of the reference administration to 1
		reg1<-lm(log(gamma)~factor(itempar)+t-1,data=tabDscrmn,x=TRUE) # estimation of the A equating coefficients with the multiple mean-geometric mean method (Haberman, 2009)
		As<-rep(1,T)
		#names(As)<-Ts
		As[-base]<-exp(reg1$coef)[(ni+1):(ni+T-1)] # store A equating coefficients
		aj<-exp(reg1$coef[1:ni]) # store synthetic discrimination parameters
		names(aj)<-substr(names(aj),16,1000)
	}
	# MULTIPLE MEAN-MEAN
	# step 1: estimation of A equating coefficients
	if (nrow(tabDscrmn)>0 & method1=="mean-mean") {
		aj1T<-reshape(tabDscrmn[,c("itempar","gamma","t")],direction="wide",v.names="gamma",timevar="t", idvar = "itempar") # discrimination parameters in wide format
		rownames(aj1T)<-aj1T$itempar
		mmout<-ipfRcpp(as.matrix(aj1T[,-1]),base,0.0000001)
		As<-mmout[[1]]
		aj_tmp<-mmout[[2]]
		#names(As)<-Ts
		names(aj_tmp)<-aj1T$itempar
		aj<-aj_tmp[names(aj)]
	}

	# if there are no discrimination parameters set all A equating coefficients to 1
	if (nrow(tabDscrmn)==0) {
		As<-rep(1,T)
		#names(As)<-Ts
		aj<-rep(1,ni)
		names(aj)<-names(seaj)<-paste("Dscrmn",substr(itms,8,100),sep=".")
	}
	
	# step 2: estimation of B equating coefficients (for both multiple mean-mean and multiple mean geometric mean methods)
	vettA<-As[tabDffclt$t]
	tabDffclt$gammaA<-tabDffclt$gamma*vettA # response variable for the second regression model
	tabDffclt$t<-as.factor(tabDffclt$t)
	tabDffclt$t <- relevel(tabDffclt$t, ref = base) #set the B equating coefficient of the reference administration to 0
	reg2<-lm(gammaA~factor(itempar)+t-1,data=tabDffclt,x=T) # estimation of the B equating coefficients
	Bs<-rep(0,T)
	#names(Bs)<-Ts
	Bs[-base]<- -reg2$coef[(ni+1):(ni+T-1)] # store B equating coefficients
	bj<-reg2$coef[1:ni] # store synthetic difficulty parameters
	names(bj)<-substr(names(bj),16,1000)

	return(c(aj,As[-base],bj,Bs[-base]))
}


# multiple equating coefficients
multiec<-function(mods, base=1, method="mean-mean", se=TRUE, nq=30, start=NULL, iter.max=100000, obsinf = TRUE)
{
  if (method=="mean-mean" | method=="mean-gmean") out<-multiec_moments(mods=mods,base=base,method=method,se=se)
  if (method=="irf" | method=="trf") out<-multiec_irf(mods=mods,base=base,method=method,se=se,nq=nq,start=start,iter.max=iter.max)
  if (method=="lik") out<-multiec_lik(mods=mods,base=base,se=se,start=start,iter.max=iter.max, obsinf = obsinf)
  return(out)
}


# multiple equating coefficients with methods based on moments
multiec_moments<-function(mods,base,method,se){
	cat("Computation of equating coefficients  .  .  .  . \n")
	itms<-c() # item labels
	for (k in 1:length(mods)) itms<-c(itms,names(mods[[k]]$coef))
	itms<-sort(unique(itms))
	itms<-itms[substr(itms,1,6)!="Gussng"]
	T<-length(mods) # number of administrations
	modsnames<-names(mods)
	tab<-c() # table with all item parameters. Columns: itempar (label of item parameters), gamma (item parameter value), t (administration)
	for (k in 1:length(mods)) {
		tab1<-data.frame(itempar=names(mods[[k]]$coef),gamma=mods[[k]]$coef,stringsAsFactors=FALSE)
		tab1$t<-k
		tab<-rbind(tab,tab1)
	}
	rownames(tab)<-paste(tab$itempar,tab$t,sep=".")
	tabDscrmn<-tab[substr(tab$itempar,1,6)=="Dscrmn",] # table with discrimination parameters
	tabDffclt<-tab[substr(tab$itempar,1,6)=="Dffclt",] # table with difficulty parameters
	ni<-length(grep("Dffclt",itms)) # numer of items
	
	# MULTIPLE MEAN-GEOMETRIC MEAN (performed anyway)
	# step 1: estimation of A equating coefficients
	if (nrow(tabDscrmn)>0) {
		if (any(tabDscrmn$gamma<0)) warning("Negative discrimination parameter. Converted to positive value.")
		tabDscrmn$gamma<-abs(tabDscrmn$gamma) # negative discriminations converted to positive values
		tabDscrmn$t<-as.factor(tabDscrmn$t)
		tabDscrmn$t <- relevel(tabDscrmn$t, ref = base) #set the A equating coefficient of the reference administration to 1
		reg1<-lm(log(gamma)~factor(itempar)+t-1,data=tabDscrmn,x=TRUE) # estimation of the A equating coefficients with the multiple mean-geometric mean method (Haberman, 2009)
		As<-rep(1,T)
		names(As)<-modsnames
		As[-base]<-exp(reg1$coef)[(ni+1):(ni+T-1)] # store A equating coefficients
		aj<-exp(reg1$coef[1:ni]) # store synthetic discrimination parameters
		names(aj)<-substr(names(aj),16,1000)
	}

	# MULTIPLE MEAN-MEAN
	# step 1: estimation of A equating coefficients
	if (nrow(tabDscrmn)>0 & method=="mean-mean") {
		aj1T<-reshape(tabDscrmn[,c("itempar","gamma","t")],direction="wide",v.names="gamma",timevar="t", idvar = "itempar") # discrimination parameters in wide format
		rownames(aj1T)<-aj1T$itempar
		mmout<-ipfRcpp(as.matrix(aj1T[,-1]),base,0.0000001)
		As<-mmout[[1]]
		aj_tmp<-mmout[[2]]
		names(As)<-modsnames
		names(aj_tmp)<-aj1T$itempar
		aj<-aj_tmp[names(aj)]
	}

	# if there are no discrimination parameters set all A equating coefficients to 1 and standard errors to 0
	if (nrow(tabDscrmn)==0) {
		As<-rep(1,T)
		names(As)<-modsnames
		seA<-rep(0,T)
		names(seA)<-modsnames
		aj<-rep(1,ni)
		seaj<-rep(0,ni)
		names(aj)<-names(seaj)<-paste("Dscrmn",substr(itms,8,100),sep=".")
	}
	
	# step 2: estimation of B equating coefficients (for both multiple mean-mean and multiple mean geometric mean methods)
	vettA<-As[tabDffclt$t]
	tabDffclt$gammaA<-tabDffclt$gamma*vettA # response variable for the second regression model
	tabDffclt$t<-as.factor(tabDffclt$t)
	tabDffclt$t <- relevel(tabDffclt$t, ref = base) #set the B equating coefficient of the reference administration to 0
	reg2<-lm(gammaA~factor(itempar)+t-1,data=tabDffclt,x=T) # estimation of the B equating coefficients
	Bs<-rep(0,T)
	names(Bs)<-modsnames
	Bs[-base]<- -reg2$coef[(ni+1):(ni+T-1)] # store B equating coefficients
	bj<-reg2$coef[1:ni] # store synthetic difficulty parameters
	names(bj)<-substr(names(bj),16,1000)

	if (se) {
		# ================================
		# computation of standard errors
		# ================================

		cat("Computation of standard errors ")
		
		tabDscrmn$t<-as.numeric(as.character(tabDscrmn$t))
		
		if (nrow(tabDscrmn)>0) {
			#VarAll<-VarExt(mods) 
			vars<-lapply(mods,FUN=function(x) x$var)
			for (i in 1:length(vars)) rownames(vars[[i]])<-colnames(vars[[i]])<-paste(rownames(vars[[i]]),i,sep=".")
			VarAll<-VarExtRcpp(vars)
			VarAllNames<-c()
			for (i in 1:length(vars)) VarAllNames<-c(VarAllNames,rownames(vars[[i]]))
			rownames(VarAll)<-colnames(VarAll)<-VarAllNames
			
			X1<-reg1$x # design matrix
			X2<-reg2$x # design matrix
			X2[,(ni+1):(ni+T-1)]<- -X2[,(ni+1):(ni+T-1)]
			colnames(X1)[1:ni]<-substr(colnames(X1)[1:ni],16,1000)
			colnames(X2)[1:ni]<-substr(colnames(X2)[1:ni],16,1000)
			sel1<-rownames(tabDscrmn)
			sel2<-rownames(tabDffclt)
			sel<-c(sel1,sel2)
			
			cat(" . ")
			
			if (method=="mean-gmean") {
				invXX1<-summary(reg1)$cov.unscaled
				#pAa<-diag(exp(reg1$coef))%*%invXX1%*%t(X1)%*%diag(1/tabDscrmn$gamma)#) # see Battauz (2016) Equation (37)
				tmp2<-matD(t(X1),1/tabDscrmn$gamma)
				tmp1<-Dmat(exp(reg1$coef),invXX1)
				pAa<-tmp1%*%tmp2
				rownames(pAa)<-rownames(t(X1))
				colnames(pAa)<-sel1
				
				cat(" . ")
				
			}
			if (method=="mean-mean") {
				pAa0<-jacobian(func=ipf4der,x=tabDscrmn$gamma,itms=tabDscrmn$itempar,t=tabDscrmn$t,base=base,aj1T=aj1T)

				cat(" . ")

				rownames(pAa0)<-modsnames
				colnames(pAa0)<-rownames(tabDscrmn)
				sumajs<-tapply(tabDscrmn$gamma,tabDscrmn$itempar,sum)
				sumAs<-tapply(As[tabDscrmn$t],tabDscrmn$itempar,sum)
				sumpAa0<-aggregate(pAa0[tabDscrmn$t,],by=list(Group.1=tabDscrmn$itempar),FUN=sum)
				rn<-sumpAa0$Group.1
				sumpAa0<-sumpAa0[,-1]
				tmp1<-matrix(NA,nrow(sumpAa0),ncol(sumpAa0))
				for (i in 1:ncol(sumpAa0)) tmp1[,i]<-sumajs/sumAs^2*sumpAa0[,i]
				rownames(tmp1)<-rn
				colnames(tmp1)<-rownames(tabDscrmn)
				tmp2<-as.matrix(1/sumAs)
				tmp2<-tmp2[,rep(1,ncol(tmp1))]
				tmp2<-tmp2*t(X1[,1:ni])
				colnames(tmp2)<-colnames(tmp1)
				pAa1<-tmp2-tmp1
				pAa<-rbind(pAa1,pAa0[-base,]) # see Battauz (2016) Equations (41) and (42)
			
			}

			cat(" . ")

			pAb<-matrix(0,ncol(X1),nrow(X1))
			colnames(pAb)<-sel2
			tmp<-Dmat(tabDffclt$gamma,-X2[,(ni+1):(ni+T-1)])
			pBa2<-t(X2)%*%tmp
			invXX2<-summary(reg2)$cov.unscaled
			colnames(invXX2)[1:ni]<-rownames(invXX2)[1:ni]<-colnames(X2)[1:ni]
			invXX2[1:ni,-(1:ni)]<- -invXX2[1:ni,-(1:ni)]
			invXX2[-(1:ni),1:ni]<- -invXX2[-(1:ni),1:ni]
			pBa<-(invXX2%*%(pBa2))%*%pAa[(ni+1):(ni+T-1),] # see Battauz (2016) Equation (39)
			#pBb_bis<-invXX2%*%t(X2)%*%diag(vettA)
			pBb<-matD(invXX2%*%t(X2),vettA) # see Battauz (2016) Equation (40)
			tmp1<-cbind(pAa,pAb)
			tmp2<-cbind(pBa,pBb)
			colnames(tmp2)<-colnames(tmp1)
			part<-rbind(tmp1,tmp2)
			cat(" . \n")
			varABgamma<-part%*%VarAll[sel,sel]%*%t(part) # covariance matrix of equating coefficients and synthetic item parameters
			seABgamma<-diag(varABgamma)^0.5
			seA<-rep(0,T)
			names(seA)<-modsnames
			seA[-base]<-seABgamma[(ni+1):(ni+T-1)] # standard errors of A equating coefficients
			seB<-rep(0,T)
			names(seB)<-modsnames
			seB[-base]<-seABgamma[(ni+ncol(X1)+1):(ni+ncol(X1)+T-1)] # standard errors of B equating coefficients
			varAB<-varABgamma[c((ni+1):(ni+T-1),(ni+ncol(X1)+1):(ni+ncol(X1)+T-1)),c((ni+1):(ni+T-1),(ni+ncol(X1)+1):(ni+ncol(X1)+T-1))] # covariance matrix of equating coefficients
			seaj<-seABgamma[1:ni] # standard errors of synthetic discrimination parameters aj
			sebj<-seABgamma[(ni+T):(2*ni+T-1)] # standard errors of synthetic difficulty parameters bj
		
			partial<-t(part[c((ni+1):(ni+T-1),(ni+ncol(X1)+1):(ni+ncol(X1)+T-1)),]) # derivatives of A and B equating coefficients with respect to the item parameters
			namesAB<-c(paste("A",(1:T)[-base],sep="."),paste("B",(1:T)[-base],sep="."))
			colnames(partial)<-namesAB
			rownames(varAB)<-colnames(varAB)<-namesAB
		}
		
		if (nrow(tabDscrmn)==0) {
			X2<-reg2$x
			mat<-summary(reg2)$cov.unscaled%*%t(X2)
			cat(" . ")
			colnames(mat)<-sel<-paste(tabDffclt$itempar,tabDffclt$t,sep=".")
			vars<-lapply(mods,FUN=function(x) x$var)
			for (i in 1:length(vars)) rownames(vars[[i]])<-colnames(vars[[i]])<-paste(rownames(vars[[i]]),i,sep=".")
			VarAll<-VarExtRcpp(vars)
			cat(" . ")
			VarAllNames<-c()
			for (i in 1:length(vars)) VarAllNames<-c(VarAllNames,rownames(vars[[i]]))
			rownames(VarAll)<-colnames(VarAll)<-VarAllNames
			sel<-paste(tabDffclt$itempar,tabDffclt$t,sep=".")
			varB<-mat%*%VarAll[sel,sel]%*%t(mat)
			cat(" . ")
			seBt<-sqrt(diag(varB))
			seB<-rep(0,T)
			names(seB)<-modsnames
			seB[-base]<-seBt[(ni+1):(ni+T-1)]
			sebj<-seBt[1:ni]
			names(sebj)<-substr(names(sebj),16,1000)
			cat(" . \n")
			partial<-matrix(0,nrow(tab),(T-1)*2)
			rownames(partial)<-paste(tab$itempar,tab$t,sep=".")
			tmp5<-t(mat[(ni+1):(ni+T-1),])
			partial[rownames(tmp5),T:(2*T-2)]<-tmp5
			namesAB<-c(paste("A",(1:T)[-base],sep="."),paste("B",(1:T)[-base],sep="."))
			colnames(partial)<-namesAB
			varAB<-matrix(0,(T-1)*2,(T-1)*2)
			varAB[T:(T*2-2),T:(T*2-2)]<-varB[(ni+1):(ni+T-1),(ni+1):(ni+T-1)]
			rownames(varAB)<-colnames(varAB)<-namesAB
		}
	}
	else {
		seA<-rep(NA,T)
		seB<-rep(NA,T)
		varAB<-matrix(NA,T,T)
		seaj<-NULL
		sebj<-NULL
		vars<-NULL
		partial<-NULL
	}
	
	tabwide<-reshape(tab[,c("itempar","gamma","t")],direction="wide",v.names="gamma",timevar="t", idvar = "itempar") # parameters in wide format
	colnames(tabwide)<-c("Item",modsnames)
	tabwide<-tabwide[order(tabwide$Item),]
	tabwides<-tabwide[,-1]
	for (t in 1:T) {
		selDffclt<-grep("Dffclt",tabwide$Item)
		selDscrmn<-grep("Dscrmn",tabwide$Item)
		tabwides[selDffclt,t]<-tabwides[selDffclt,t]*As[t]+Bs[t] # conversion of difficulty parameters to the scale of the base form
		tabwides[selDscrmn,t]<-tabwides[selDscrmn,t]/As[t] # conversion of discrimination parameters to the scale of the base form
	}
	colnames(tabwides)<-paste(colnames(tabwides),modsnames[base],sep=".as.")
	tabwides<-tabwides[,-base]
	tabwide<-cbind(tabwide,tabwides)
	itmp<-2
	if (sum(substr(tabwide$Item,1,6)=="Dscrmn")==0) itmp=1
	if (sum(substr(tabwide$Item,1,6)=="Gussng")>0) itmp=3

	out<-list(A=As,B=Bs,se.A=seA,se.B=seB,varAB=varAB,as=aj,bs=bj,se.as=seaj,se.bs=sebj,tab=tabwide,varFull=vars,partial=partial,itmp=itmp,method=method,basename=modsnames[base])
	class(out)<-"mlteqc"
	return(out)
}


print.mlteqc<-function(x, ...)
{
	cat("Multiple equating coefficients \n")
	cat("Method: ")
	cat(x$method,"\n")
}


summary.mlteqc <- function(object, ...)
{
	ct<-data.frame(EQ=rep(c("A","B"),each=length(object$A)),Form=c(names(object$A),names(object$B)),Estimate=c(object$A,object$B),StdErr=c(object$se.A,object$se.B))
	cat("Equating coefficients:\n")
	print(ct,digits=5,row.names=F)
}




# product of a matrix mat and a diagonal matrix with diagonal D: mat%*%diag(D)
matD<-function(mat,D){
	nr<-nrow(mat)
	for (i in 1:nr) mat[i,]<-mat[i,]*D
	return(mat)
	rownames(mat)<-colnames(mat)<-NULL
}	
		
# product of a diagonal matrix with diagonal D and a matrix mat: diag(D)%*%mat
Dmat<-function(D,mat){
	nc<-ncol(mat)
	for (i in 1:nc) mat[,i]<-D*mat[,i]
	rownames(mat)<-colnames(mat)<-NULL
	return(mat)
}			





itm.mlteqc<-function(x, ...) x$tab


eqc.mlteqc<-function(x, ...)
{
	A1<-x$A
	B1<-x$B
	out<-data.frame(A=A1,B=B1)
	rownames(out)<-NULL
	return(out)
}



score.mlteqc<-function(obj, method="TSE", D=1, scores=NULL, se=TRUE, nq=30, w=0.5, theta=NULL, weights=NULL, ...)
{
  if (!is.null(scores)) if (any(round(scores)!=scores)) stop("Scores should be integer values.")
  
  itmpar<-itm(obj)
  basename<-obj$basename
  modsnames<-names(obj$A)
  T<-length(obj$A)
  base<-(1:T)[modsnames==basename]
  out<-NULL
  for (t in 1:T) {
    if (base!=t) {
      sel<-c("Item",modsnames[t],basename,paste(modsnames[t],basename,sep=".as."))
      sel1<-c(paste("A",t,sep="."),paste("B",t,sep="."))
      itm_prepare<-equateIRT:::score_prepare(itmpar[,sel],suff1=paste(".",t,sep=""),suff2=paste(".",base,sep=""))
      out_t<-equateIRT:::score_compute(method=method,itm_prepare=itm_prepare,D=D,varFull=obj$varFull,partial=obj$partial[,sel1],varAB=obj$varAB[sel1,sel1],itmp=obj$itmp,A=obj$A[t],B=obj$B[t],scores=scores,se=se,nq=nq,w=w,theta=theta,weights=weights,names=sel[3:4])
      if (is.null(out)) {
				k<-ncol(out_t)
				if (se) colnames(out_t)[k]<-paste(colnames(out_t)[k],colnames(out_t)[k-1],sep="_")
				out<-out_t
			}
      else {
        if (method=="TSE") out_t<-subset(out_t,select=3:ncol(out_t))
        if (method=="OSE") out_t<-subset(out_t,select=2:ncol(out_t))
        if (se) colnames(out_t)[2]<-paste(colnames(out_t)[2],colnames(out_t)[1],sep="_")
        out<-cbind(out,out_t)
      }
    }
  }
  return(out)
}






# ProfLik<-function(par,itmpar,itmvar,num.forms,base,DffcltNum,DscrmnNum)
# {
#   coefs<-NULL
#   A<-rep(1,num.forms)
#   A[-base]<-par[1:(num.forms-1)]
#   B<-rep(0,num.forms)
#   B[-base]<-par[((num.forms-1)+1):(2*(num.forms-1))]
#   itmpar$A<-A[itmpar$t]
#   itmpar$B<-B[itmpar$t]
#   itmpar$Y<-9999
#   itmpar[DscrmnNum,"Y"]<-itmpar[DscrmnNum,coef/A]
#   itmpar[DffcltNum,"Y"]<-itmpar[DffcltNum,coef*A+B]
#   Y<-itmpar$Y
#   X<-model.matrix(~factor(items)-1,data=itmpar)
#   row.names(X)<-itmpar$items.t
#   X_list<-list()
#   for (f in 1:num.forms) X_list[[f]]<-X[itmpar$t==f,]
#   Y_list<-list()
#   for (f in 1:num.forms) Y_list[[f]]<-Y[itmpar$t==f]
#   omega<-itmvar
#   for (t in 1:num.forms)
#   {
#     selds<-grep("Dscrmn",rownames(omega[[t]]))
#     seldf<-grep("Dffclt",rownames(omega[[t]]))
#     omega[[t]][selds,selds]<-omega[[t]][selds,selds]/A[t]^2
#     omega[[t]][seldf,seldf]<-omega[[t]][seldf,seldf]*A[t]^2
#   }
#   V.chol <- try(lapply(omega, Matrix::chol))
#   V.chol.inv <- try(lapply(V.chol, Matrix::solve))
#   if (!isa(V.chol.inv,"try-error"))
#   {
#     t.V.chol.inv.X<-mapply(FUN=Matrix::crossprod,V.chol.inv,X_list,SIMPLIFY=F)
#     t.V.chol.inv.X<-do.call("rbind", t.V.chol.inv.X)
#     tX.V.chol.inv_t.V.chol.inv.X<-Matrix::crossprod(t.V.chol.inv.X)
#     # max(abs(tX.V.chol.inv_t.V.chol.inv.X-tX.Oinv.X))
#     t.V.chol.inv.Y<-mapply(FUN=Matrix::crossprod,V.chol.inv,Y_list,SIMPLIFY=F)
#     t.V.chol.inv.Y<-do.call("rbind", t.V.chol.inv.Y)
#     tX.V.chol.inv_t.V.chol.inv.Y<-Matrix::crossprod(t.V.chol.inv.X,t.V.chol.inv.Y)
#     beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X,tX.V.chol.inv_t.V.chol.inv.Y)
#     # max(abs(beta-beta1))
#     
#     rownames(beta)<-substr(colnames(X),14,100)
#     itmpar$coefs<-beta[itmpar$items,]
#     itmpar$mean<-9999
#     itmpar[DscrmnNum,mean:=coefs*A]
#     itmpar[DffcltNum,mean:=(coefs-B)/A]
#     logden<-itmpar[,dmvnorm(x=coef,mean=mean,sigma=itmvar[[t]],log=TRUE),by=t]
#     out<-sum(logden$V1)
#   }
#   else out<- 1e100000000
#   -out
# }
# 
# 
# 
# ProfLik_1PL<-function(par,itmpar,itmvar,num.forms,base)
# {
#   coefs<-NULL
#   B<-rep(0,num.forms)
#   B[-base]<-par
#   itmpar$B<-B[itmpar$t]
#   itmpar$Y<-9999
#   itmpar$Y<-itmpar[,coef+B]
#   Y<-itmpar$Y
#   X<-model.matrix(~factor(items)-1,data=itmpar)
#   X_list<-list()
#   for (f in 1:num.forms) X_list[[f]]<-X[itmpar$t==f,]
#   Y_list<-list()
#   for (f in 1:num.forms) Y_list[[f]]<-Y[itmpar$t==f]
#   omega<-itmvar
#   V.chol <- try(lapply(omega, Matrix::chol))
#   V.chol.inv <- try(lapply(V.chol, Matrix::solve))
#   if (!isa(V.chol.inv,"try-error"))
#   {
#     t.V.chol.inv.X<-mapply(FUN=Matrix::crossprod,V.chol.inv,X_list,SIMPLIFY=F)
#     t.V.chol.inv.X<-do.call("rbind", t.V.chol.inv.X)
#     tX.V.chol.inv_t.V.chol.inv.X<-Matrix::crossprod(t.V.chol.inv.X)
#     # max(abs(tX.V.chol.inv_t.V.chol.inv.X-tX.Oinv.X))
#     t.V.chol.inv.Y<-mapply(FUN=Matrix::crossprod,V.chol.inv,Y_list,SIMPLIFY=F)
#     t.V.chol.inv.Y<-do.call("rbind", t.V.chol.inv.Y)
#     tX.V.chol.inv_t.V.chol.inv.Y<-Matrix::crossprod(t.V.chol.inv.X,t.V.chol.inv.Y)
#     beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X,tX.V.chol.inv_t.V.chol.inv.Y)
#     # max(abs(beta-beta1))
#     
#     rownames(beta)<-substr(colnames(X),14,100)
#     itmpar$coefs<-beta[itmpar$items,]
#     itmpar$mean<-9999
#     itmpar[,mean:=(coefs-B)]
#     logden<-itmpar[,dmvnorm(x=coef,mean=mean,sigma=itmvar[[t]],log=TRUE),by=t]
#     out<-sum(logden$V1)
#   }
#   else out<- 1e100000000
#   -out
# }
# 



# estimates of the item parameters on a common scale using the estimates from all the forms
getBetas<-function(par,itmpar,itmvar,num.forms,base,DffcltNum,DscrmnNum,X_list)
{
  A<-rep(1,num.forms)
  A[-base]<-par[1:(num.forms-1)]
  B<-rep(0,num.forms)
  B[-base]<-par[((num.forms-1)+1):(2*(num.forms-1))]
  itmpar$A<-A[itmpar$t]
  itmpar$B<-B[itmpar$t]
  itmpar$Y<-9999
  itmpar[DscrmnNum,"Y"]<-itmpar[DscrmnNum,coef/A]
  itmpar[DffcltNum,"Y"]<-itmpar[DffcltNum,coef*A+B]
  Y<-itmpar$Y
  Y_list<-list()
  for (f in 1:num.forms) Y_list[[f]]<-Y[itmpar$t==f]
  
  omega<-itmvar
  for (t in 1:num.forms)
  {
    selds<-grep("Dscrmn",rownames(omega[[t]]))
    seldf<-grep("Dffclt",rownames(omega[[t]]))
    omega[[t]][selds,selds]<-omega[[t]][selds,selds]/A[t]^2
    omega[[t]][seldf,seldf]<-omega[[t]][seldf,seldf]*A[t]^2
  }
  V.chol <- try(lapply(omega, Matrix::chol))
  V.chol.inv <- try(lapply(V.chol, Matrix::solve))
  V.inv <- lapply(V.chol.inv, Matrix::tcrossprod)

  if (!isa(V.chol.inv,"try-error"))
  {
    t.V.chol.inv.X<-mapply(FUN=Matrix::crossprod,V.chol.inv,X_list,SIMPLIFY=F)
    t.V.chol.inv.X<-do.call("rbind", t.V.chol.inv.X)
    
    # the following is necessary to compute standard errors of as and bs
    V.inv.X<-mapply(FUN=Matrix::crossprod,V.inv,X_list,SIMPLIFY=F)
    V.inv.X<-do.call("rbind", V.inv.X)
    
    tX.V.chol.inv_t.V.chol.inv.X<-Matrix::crossprod(t.V.chol.inv.X)
    t.V.chol.inv.Y<-mapply(FUN=Matrix::crossprod,V.chol.inv,Y_list,SIMPLIFY=F)
    t.V.chol.inv.Y<-do.call("rbind", t.V.chol.inv.Y)
    tX.V.chol.inv_t.V.chol.inv.Y<-Matrix::crossprod(t.V.chol.inv.X,t.V.chol.inv.Y)
    beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X, tX.V.chol.inv_t.V.chol.inv.Y)
  }
  else beta<-NA
  list(beta=beta, tX.Vinv.X = tX.V.chol.inv_t.V.chol.inv.X, V.inv.X = V.inv.X)
}



getBetas_1PL<-function(par,itmpar,itmvar,num.forms,base,X_list)
{
  B<-rep(0,num.forms)
  B[-base]<-par
  itmpar$B<-B[itmpar$t]
  itmpar$Y<-9999
  itmpar$Y<-itmpar[,coef+B]
  Y<-itmpar$Y
  Y_list<-list()
  for (f in 1:num.forms) Y_list[[f]]<-Y[itmpar$t==f]
  
  omega<-itmvar
  V.chol <- try(lapply(omega, Matrix::chol))
  V.chol.inv <- try(lapply(V.chol, Matrix::solve))
  V.inv <- lapply(V.chol.inv, Matrix::tcrossprod)
  
  if (!isa(V.chol.inv,"try-error"))
  {
    t.V.chol.inv.X<-mapply(FUN=Matrix::crossprod,V.chol.inv,X_list,SIMPLIFY=F)
    t.V.chol.inv.X<-do.call("rbind", t.V.chol.inv.X)
    
    # the following is necessary to compute standard errors of as and bs
    V.inv.X<-mapply(FUN=Matrix::crossprod,V.inv,X_list,SIMPLIFY=F)
    V.inv.X<-do.call("rbind", V.inv.X)
    
    tX.V.chol.inv_t.V.chol.inv.X<-Matrix::crossprod(t.V.chol.inv.X)
    t.V.chol.inv.Y<-mapply(FUN=Matrix::crossprod,V.chol.inv,Y_list,SIMPLIFY=F)
    t.V.chol.inv.Y<-do.call("rbind", t.V.chol.inv.Y)
    tX.V.chol.inv_t.V.chol.inv.Y<-Matrix::crossprod(t.V.chol.inv.X,t.V.chol.inv.Y)
    beta <- Matrix::solve(tX.V.chol.inv_t.V.chol.inv.X,tX.V.chol.inv_t.V.chol.inv.Y)
  }
  else beta<-NA
  list(beta=beta, tX.Vinv.X = tX.V.chol.inv_t.V.chol.inv.X, V.inv.X = V.inv.X)
  
}





getitmpar<-function(mods,t)
{
  ct<-mods$coefficients
  itemstype<-substr(names(ct),1,6)
  noGussng<-itemstype!="Gussng"
  ct<-ct[noGussng]
  itemstype<-itemstype[noGussng]
  itemslab<-substr(names(ct),8,100)
  data.frame(items=names(ct),itemstype=itemstype,itemslab=itemslab,coef=ct,t=t)
}

delGussng<-function(x)
{
  which_guess<-grep("Gussng",rownames(x))
  if (length(which_guess)>0) x<-x[-which_guess,-which_guess]
  x
}

multiec_lik <- function(mods, base, se, obsinf, start, iter.max)
{
  itemstype<-items<-Y<-items.t<-NULL
  cat("Computation of equating coefficients  .  .  . ")
  num.forms<-length(mods)
  
  itmpar<-mapply(FUN=getitmpar,mods=mods,t=1:num.forms,SIMPLIFY=FALSE)
  itmpar<-rbindlist(itmpar)
  itmpar[,items.t:=paste(items,t,sep=".")]
  
  DscrmnNum<-grep("Dscrmn",itmpar$items)
  DffcltNum<-grep("Dffclt",itmpar$items)
  
  itmvar<-lapply(mods,FUN=function(x) x$var)
  itmvar<-lapply(itmvar,delGussng)
  for (i in 1:length(itmvar)) rownames(itmvar[[i]])<-colnames(itmvar[[i]])<-paste(rownames(itmvar[[i]]),i,sep=".")

  itmp<-2
  if (sum(substr(itmpar$items,1,6)=="Dscrmn")==0) itmp=1
  if (sum(substr(itmpar$itmes,1,6)=="Gussng")>0) itmp=3
  
  if (itmp>=2)
  {
    if (inherits(start,"mlteqc")) start<-c(start$A[-base],start$B[-base])
    if (is.null(start))
    {
      capture.output(irf<-multiec(mods,se=FALSE,method = "irf",base=base))
      start<- c(irf$A[-base],irf$B[-base])
    }
  }
  if (itmp==1)
  {
    if (inherits(start,"mlteqc")) start<-start$B[-base]
    if (is.null(start))
    {
      irf<-multiec(mods,se=FALSE,method = "irf",base=base)
      start<- irf$B[-base]
    }
  }
  
  ini<-start
  
  se.A<-se.B<-rep(NA,num.forms)
  varAB<-NULL

  X<-model.matrix(~factor(items)-1,data=itmpar)
  row.names(X)<-itmpar$items.t
  colnames(X)<-substr(colnames(X),14,100)
  X_list<-list()
  for (f in 1:num.forms) X_list[[f]]<-X[itmpar$t==f,]
  pos<-match(itmpar$items,colnames(X))
  
  if (itmp>=2) opt<-try(optim(par=ini,fn=profLikRcpp,
                              coef=itmpar$coef, t=itmpar$t-1, X_list=X_list,
                              itmvar=itmvar, numforms=num.forms, notbase=(1:num.forms)[-base]-1,
                              DffcltNum=DffcltNum-1,DscrmnNum=DscrmnNum-1,pos=pos-1,
                              method="BFGS",hessian=se,control=list(maxit=iter.max)))
  if (itmp==1) opt<-try(optim(par=ini,fn=profLikRcpp_1PL,
                              coef=itmpar$coef, t=itmpar$t-1, X_list=X_list,
                              itmvar=itmvar, numforms=num.forms, 
                              notbase=(1:num.forms)[-base]-1,pos=pos-1,
                              method="BFGS",hessian=se,control=list(maxit=iter.max)))
  cat(" .   \n")
  if (!isa(opt,"try-error"))
  {
    A<-rep(1,num.forms)
    B<-rep(0,num.forms)
    if (itmp>=2) 
    {
      A[-base]<-opt$par[1:(num.forms-1)]
      B[-base]<-opt$par[num.forms:(2*(num.forms-1))]
    }
    if (itmp==1)
    {
      B[-base]<-opt$par
    }
    modsnames<-names(mods)
    names(A)<-names(B)<-modsnames
    conv<-opt$convergence
    
    if (itmp>=2)
    {
      outBetas<-getBetas(par=opt$par,itmpar=itmpar,itmvar=itmvar,num.forms=num.forms,base=base,DffcltNum,DscrmnNum,X_list=X_list)
      betas<-outBetas$beta
      as<-betas[grep("Dscrmn",rownames(betas)),]
      bs<-betas[grep("Dffclt",rownames(betas)),]
    }
    if (itmp==1)
    {
      outBetas<-getBetas_1PL(par=opt$par,itmpar=itmpar,itmvar=itmvar,num.forms=num.forms,base=base,X_list=X_list)
      betas<-outBetas$beta
      rownames(betas)<-substr(colnames(X),14,100)
      bs<-betas[grep("Dffclt",rownames(betas)),]
      as<-NULL
    }
    
    itmpar$A<-A[itmpar$t]
    itmpar$B<-B[itmpar$t]
    itmpar$Y<-9999
    itmpar[itemstype=="Dscrmn",Y:=coef/A]
    itmpar[itemstype=="Dffclt",Y:=coef*A+B]
    
    tab<-data.table::dcast(itmpar, items ~t,value.var = c("coef","Y"))
    tab<-as.data.frame(tab)
    sel<-which(colnames(tab)==paste("Y",base,sep="_"))
    tab<-tab[,-sel] # delete base form converted
    colnames(tab)<-c("Item",modsnames,paste(modsnames[-base],modsnames[base],sep=".as."))
    
    partial<-NULL

    if (se)
    {
      cat("Computation of standard errors ")
      if (obsinf)
      {
        varAB<-Matrix::solve(opt$hessian)
        cat(" .  .  .  . \n")
      }
      else
      {
        cat(" . ")
        derS_AB<-opt$hessian
        common<-as.character(itmpar$items[duplicated(itmpar$items)])
        itmpar_common<-itmpar[items%in%common]
        gamma<-itmpar_common$coef
        names(gamma)<-itmpar_common$items.t
        if (itmp>=2) derS_gamma<-jacobian(func=derAB,x=gamma,par=opt$par,itmpar=itmpar,itmvar=itmvar,num.forms=num.forms,method="simple",base=base,DffcltNum=DffcltNum,DscrmnNum=DscrmnNum,X_list=X_list,pos=pos)
        if (itmp==1) derS_gamma<-jacobian(func=derAB_1PL,x=gamma,par=opt$par,itmpar=itmpar,itmvar=itmvar,num.forms=num.forms,method="simple",base=base,X_list=X_list,pos=pos)
        cat(" . ")
        colnames(derS_gamma)<-itmpar_common$items.t
        derAB_gamma <- -Matrix::solve(derS_AB)%*%derS_gamma
        
        partial<-t(derAB_gamma) # derivatives of A and B equating coefficients with respect to the item parameters
        if (itmp==1) partial<-cbind(matrix(0,nrow(partial),num.forms-1),partial) 
        namesAB<-c(paste("A",(1:num.forms)[-base],sep="."),paste("B",(1:num.forms)[-base],sep="."))
        colnames(partial)<-namesAB

        var_gamma<-bdiag(itmvar)
        
        VarNames<-c()
        for (i in 1:length(itmvar)) VarNames<-c(VarNames,rownames(itmvar[[i]]))
        rownames(var_gamma)<-colnames(var_gamma)<-VarNames

        sel<-colnames(derAB_gamma)
        cat(" . ")
        varAB<-derAB_gamma %*% var_gamma[sel,sel] %*% t(derAB_gamma)
        cat(" . \n")
      }
      se.A<-se.B<-rep(0,num.forms)
      names(se.A)<-names(se.B)<-modsnames
      if (itmp>=2)
      {
        se.A[-base]<-Matrix::diag(varAB)[1:(num.forms-1)]^0.5
        se.B[-base]<-Matrix::diag(varAB)[(num.forms):(2*num.forms-2)]^0.5
        namesAB<-c(paste("A",(1:num.forms)[-base],sep="."),paste("B",(1:num.forms)[-base],sep="."))
        colnames(varAB)<-rownames(varAB)<-namesAB
      }
      if (itmp==1)
      {
        se.B[-base]<-Matrix::diag(varAB)^0.5
        tmp<-matrix(0,(2*num.forms-2),(2*num.forms-2))
        tmp[num.forms:(2*num.forms-2),num.forms:(2*num.forms-2)]<-as.matrix(varAB)
        varAB<-tmp
        namesAB<-c(paste("A",(1:num.forms)[-base],sep="."),paste("B",(1:num.forms)[-base],sep="."))
        colnames(varAB)<-rownames(varAB)<-namesAB
      }
      
      # computation of se of as and bs
      tX.Vinv.X<-outBetas$tX.Vinv.X
      V.inv.X<-outBetas$V.inv.X
      der_beta_gamma<-Matrix::tcrossprod(solve(tX.Vinv.X),V.inv.X)
      
      if (obsinf)
      {
        var_gamma<-bdiag(itmvar)
        VarNames<-c()
        for (i in 1:length(itmvar)) VarNames<-c(VarNames,rownames(itmvar[[i]]))
        rownames(var_gamma)<-colnames(var_gamma)<-VarNames
      }
      sel<-colnames(der_beta_gamma)
      var_ab<-der_beta_gamma %*% var_gamma[sel,sel] %*% t(der_beta_gamma)
      se_ab<-Matrix::diag(var_ab)^0.5
      se.as<- se_ab[grep("Dscrmn",names(se_ab))]
      se.bs<- se_ab[grep("Dffclt",names(se_ab))]
    }
  }
  else
  {
    A<-B<-rep(NA,length(se.A))
    conv<-"optimization failed"
    as<-bs<-NULL
  }
  out<-list(A=A,B=B,se.A=se.A,se.B=se.B,varAB=varAB,as=as,bs=bs,se.as=se.as,
            se.bs=se.bs,tab=tab,varFull=itmvar,partial=partial,itmp=itmp,method="lik",
            basename=modsnames[base],convergence=conv)
  class(out)<-"mlteqc"
  return(out)
}

# this function is used for finding the second derivatives of the profile log-likelihood
# with respect to the equating coefficients and the item parameters
# x = item parameters
derAB<-function(x,par,itmpar,itmvar,num.forms,base,DffcltNum,DscrmnNum,X_list,pos)
{
  itmpar$coef<-x[itmpar$items.t]
  grad(func=profLikRcpp,x=par,coef=itmpar$coef, t=itmpar$t-1, X_list=X_list,itmvar=itmvar, numforms=num.forms, notbase=(1:num.forms)[-base]-1,DffcltNum=DffcltNum-1,DscrmnNum=DscrmnNum-1,pos=pos-1)
}

derAB_1PL<-function(x,par,itmpar,itmvar,num.forms,base,X_list,pos)
{
  itmpar$coef<-x[itmpar$items.t]
  grad(func=profLikRcpp_1PL,x=par,coef=itmpar$coef, t=itmpar$t-1, X_list=X_list,itmvar=itmvar, numforms=num.forms,notbase=(1:num.forms)[-base]-1,pos=pos-1)
}


summary.likeqc<-function (object, ...)
{
  ct<-data.frame(EQ=rep(c("A","B"),each=length(object$A)),Form=c(names(object$A),names(object$B)),Estimate=c(object$A,object$B),StdErr=c(object$se.A,object$se.B))
  cat("Equating coefficients:\n")
  print(ct,digits=5,row.names=F)
}





