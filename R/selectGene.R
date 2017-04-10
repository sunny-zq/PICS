#' @title Gene selection using SPLS model
#'
#' @description Sparse Partial Least Square (SPLS) step for gene selection and dimension reduction.By applying SPLS to each pathway, we achieve the goal of gene selection and dimension reduction at the same time.
#'
#' @param object output list of prefilter step
#' @param fold The number of folds to use to perform the cross-validation process.
#' @param K the maximum number of hidden features in spls.
#' @param etas Thresholding parameter. eta should be between 0 and 1.
#' @param seed random seed that was set,  default = 123.
#'
#' @import stats foreach plsRcox survival methods
#' @export
#' @docType methods
#' @rdname selectGene-methods
#' @aliases selectGene
#' @aliases selectGene,Prefiltered-method
#'
#' @examples
#' data(TCGA)
#' prefilter.results=prefilter( data=TCGA$geneexpr, time=TCGA$t, status=TCGA$d, plist=TCGA$pathList )
#' gene.results=selectGene( prefilter.results )

#nt:maximum number of latent variables to consider
#se1: whether or not to use the 1se criteria for cv
#method: auc or partial likelihood deviance
#par:whether or not to do parallel computing
setMethod(
  f="selectGene",
  signature="Prefiltered",
  definition=function( object, nt=5, etas=seq(0.1,0.9,0.1), fold=5, 
  se1=TRUE, method="auc", par=FALSE, foldid=NULL, seed=123 ) {

    t<-object@inputdata$time
    d<-object@inputdata$status
    n<-length(t)
    data<-object@xlist
    pathways<-object@inputdata$pathway
    dimx=unlist( lapply(data,function(x){ncol(as.matrix(x))}) )
    
    if(is.null(foldid)){
    	set.seed(seed)
		foldid = sample(rep(seq(1:fold), length=n))
	}
	
    k.opt<-eta.opt<-NULL
    score<-genes<-beta<-w<-list()

    for(j in 1:length(pathways)){
    	
      xx<-as.matrix( data[[j]],nrow=n,ncol=dimx[j] )
      kmax<-min( nt, ncol(xx) )

	  cvs=cv_splscox(x=xx,t=t,d=d,foldid=foldid,K=kmax,eta.vec=etas,method=method,parallel=par)
  	  
  	  if(se1==TRUE){
  	  	k.opt[j]=cvs$opt.k[2]
      	eta.opt[j]=cvs$opt.eta[2] 
      }	
      ##if auc is used, this is the maximum auc criteria
      ##if partial likelihood is used, this is the min deviance criteria
  	  if(se1==FALSE){
  	  	k.opt[j]=cvs$opt.k[1]
      	eta.opt[j]=cvs$opt.eta[1] 
      }	
      
  	  cox<-coxph(Surv(t,d)~1)
      devres<-residuals(cox,type="deviance")

      spls.mod<-spls.cox (x=xx, y=devres, K=k.opt[j], eta=eta.opt[j],
                           kappa=0.5, select="pls2", scale.x=T, scale.y=F, trace=F)

      score[[j]]<- data.frame(spls.mod$plsmod$variates$X)
      beta[[j]]<-data.frame( colnames(xx),spls.mod$betahat )
      rownames(beta[[j]])<-NULL

      ##objects saved for prediction
      genes[[j]]<-colnames(xx)[spls.mod$A]
      w[[j]]<-spls.mod$pred$w

	  print(sprintf("Interation: %d", j)) 
	  flush.console()      
    }

    names(genes)<-pathways
    names(beta)<-pathways

    methods::new( "FitGene",
      score=score,
      geneSelected=genes,
      fit = list( coef=beta, direction=w ),
      dataPrefiltered=data,
      inputdata = object@inputdata
    )
  }
)
