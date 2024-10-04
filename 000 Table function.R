
library(magrittr)
library(RCurl)
library(XML)
library(survival)
library(rtf)

.ROUND <- function(X,digits=0) {
  x=rep(NA,length(X))
  for (a in 1:length(X)) {
    if (is.na(X[a])) {x[a]=NA} else {
      x[a]=formatC(X[a],format="f",digits=digits)
    }
  }
  return(x)
}

.ROUND.p <- function(X,p.digits=3) {
  p=rep(NA,length(X))
  for (a in 1:length(X)) {
    if (!is.na(X[a])) {
      if (p.digits==0) {p[a]="<1"} else {
        p[a]=as.character(round(X[a],p.digits))
        if (p[a]=="0") {
          if (p.digits==1) {p[a]="<0.1"} else {
            p[a]="<0."
            for (l in 1:(p.digits-1)) {
              p[a]=paste0(p[a],0)
            }
            p[a]=paste0(p[a],1)  
          }
        } else {
          p[a]=.ROUND(as.numeric(X[a]),p.digits)
        }
      }
    }
  }
  return(p)
}

Table1 <- function(X,Y.matrix,digits=1,digits.per=1,p.digits=3,Factor=NULL,x.name="Group",SD=TRUE,nonparametric=TRUE,Transpose=FALSE) {
  if (is.null(Factor)) {
    Factor=numeric()
    for (k in 1:ncol(Y.matrix)) {
      if (is.factor(Y.matrix[,k])) {
        Factor=c(Factor,1)
      } else {
        Factor=c(Factor,0)
      }   
    }
  }
  n.row=numeric()
  row.names=NULL
  for (k in 1:ncol(Y.matrix)) {
    if (Factor[k]==1) {
      n.row=c(n.row,1+length(levels(Y.matrix[,k])))
      row.names=c(row.names,colnames(Y.matrix)[k],paste0(colnames(Y.matrix)[k],":",levels(Y.matrix[,k])))
    } else {
      n.row=c(n.row,1)
      row.names=c(row.names,colnames(Y.matrix)[k])
    }   
  }
  Table=matrix("",nrow=sum(n.row),ncol=length(levels(factor(X)))+1)
  colnames(Table)=c(paste0(x.name,":",levels(factor(X))),"p-value")
  rownames(Table)=row.names
  position=cumsum(n.row)
  #######################################
  for (k in 1:ncol(Y.matrix)) {
    if (Factor[k]==0) {
      n.sample=numeric()     
      for (i in 1:length(levels(factor(X)))) {
        n.sample=c(n.sample,sum(is.na(Y.matrix[factor(X)==levels(factor(X))[i],k])==F))
        m=mean(Y.matrix[factor(X)==levels(factor(X))[i],k],na.rm=T)
        if (SD==TRUE) {
          s=sd(Y.matrix[factor(X)==levels(factor(X))[i],k],na.rm=T)
        } else {
          s=sd(Y.matrix[factor(X)==levels(factor(X))[i],k],na.rm=T)/sqrt(sum(!is.na(Y.matrix[factor(X)==levels(factor(X))[i],k])))
        }
        m=.ROUND(m,digits)
        s=.ROUND(s,digits)
        Table[position[k],i]=paste0(m,"±",s)
      }
      if (length(levels(factor(X)))>1) {
        p <- ''
        if (min(n.sample) > 0) {
          if (nonparametric==TRUE) {
            if (min(n.sample)>=25) {
              p=.ROUND.p(anova(lm(Y.matrix[,k]~factor(X)))$"Pr(>F)"[1],p.digits)
            } else {
              if (length(levels(factor(X)))>=3) {
                p=.ROUND.p(kruskal.test(Y.matrix[,k],factor(X))$p.value,p.digits)
              } else {
                p=.ROUND.p(wilcox.test(Y.matrix[,k]~factor(X),correct=FALSE)$p.value,p.digits)
              }
            }
            if (min(n.sample)<25) {p=paste0(p,"#")} 
          } else {
            p=.ROUND.p(anova(lm(Y.matrix[,k]~factor(X)))$"Pr(>F)"[1],p.digits)
          }
        }
        Table[position[k],length(levels(factor(X)))+1] = p
      }
    } else {
      for (i in 1:length(levels(factor(X)))) {
        for (j in 1:length(levels(Y.matrix[,k]))) {
          n=sum(Y.matrix[factor(X)==levels(factor(X))[i],k]==levels(Y.matrix[,k])[j],na.rm=T)
          prop=n/length(Y.matrix[factor(X)==levels(factor(X))[i]&is.na(Y.matrix[,k])==F&is.na(X)==F,k])
          prop=.ROUND(prop*100,digits.per)
          prop=paste0("(",prop,"%)")
          Table[position[k]-length(levels(Y.matrix[,k]))+j,i]=paste0(n,prop)
        }
      }
      if (length(levels(factor(X)))>1) {
        if (nonparametric==TRUE) {
          logic=tryCatch(chisq.test(table(Y.matrix[,k],factor(X))),error=function(e) e, warning=function(w) w)
          sig.test=try(fisher.test(table(Y.matrix[,k],factor(X))),silent=T)
          if (is(logic,"warning")&is(sig.test)!="try-error") {
            p=.ROUND.p(sig.test$p.value,p.digits)
          } else {
            p=.ROUND.p(chisq.test(table(Y.matrix[,k],factor(X)),correct=FALSE)$p.value,p.digits)
          }
          if (is(logic,"warning")&is(sig.test)!="try-error") {p=paste0(p,"#")}
        } else {
          p=.ROUND.p(chisq.test(table(Y.matrix[,k],factor(X)),correct=FALSE)$p.value,p.digits)
        }
        Table[position[k]-length(levels(Y.matrix[,k])),length(levels(factor(X)))+1]=p
      }
    }
  }
  Table.Final=cbind(rownames(Table),Table)
  colnames(Table.Final)[1]="Variable"
  rownames(Table.Final)=1:nrow(Table)
  if (length(levels(factor(X)))>1) {return(Table.Final)} else {
    Table.Final[,3]=""
    colnames(Table.Final)[2]="Total"
    return(Table.Final[,-3])
  }
}

Table2 <- function(X.matrix,Y,Adj.matrix=NULL,time=NULL,type,sig=0.05,digits=2,p.digits=3,Factor=NULL,del=" to ") {
  #if (length(X.matrix)==length(Y)) {X.matrix=matrix(X.matrix,nrow=length(X.matrix),ncol=1)}
  if ((class(X.matrix)!="matrix"|class(X.matrix)!="data.frame")&length(X.matrix)==length(Y)) {X.matrix=matrix(X.matrix,nrow=length(X.matrix),ncol=1)}
  if (length(Adj.matrix)==length(Y)) {Adj.matrix=matrix(Adj.matrix,nrow=length(Adj.matrix),ncol=1)}  
  if (is.null(Factor)) {
    Factor=numeric()
    for (k in 1:ncol(X.matrix)) {
      if (is.factor(X.matrix[,k])) {
        Factor=c(Factor,1)
      } else {
        Factor=c(Factor,0)
      }   
    }
  }
  n.row=numeric()
  row.names=NULL
  for (k in 1:ncol(X.matrix)) {
    if (Factor[k]==1) {
      n.row=c(n.row,1+length(levels(X.matrix[,k])))
      row.names=c(row.names,colnames(X.matrix)[k],paste0(colnames(X.matrix)[k],":",levels(X.matrix[,k])))
    } else {
      n.row=c(n.row,1)
      row.names=c(row.names,colnames(X.matrix)[k])
    }   
  }
  
  if (is.null(Adj.matrix)) {
    Table=matrix("",nrow=sum(n.row),ncol=3)
    if (type=="linear") {coef.name="beta"}
    if (type=="logistic") {coef.name="OR"}
    if (type=="survival") {coef.name="HR"}  
    colnames(Table)=c(paste0(coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","n")
    rownames(Table)=row.names
    position=cumsum(n.row)
  } else {
    Table=matrix("",nrow=sum(n.row),ncol=8)
    if (type=="linear") {coef.name="beta"}
    if (type=="logistic") {coef.name="OR"}
    if (type=="survival") {coef.name="HR"}  
    colnames(Table)=c(paste0("Subset 1 Crude-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","n (Subset 1)",paste0("Subset 2 Crude-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value",paste0("Subset 2 Adj-",coef.name," (",round((1-sig)*100,2),"% CI)#"),"p-value","n (Subset 2)")
    rownames(Table)=row.names
    position=cumsum(n.row)  
  }
  
  #######################################
  for (k in 1:ncol(X.matrix)) {
    if (Factor[k]==0) {
      if (type=="linear") {
        n.1=sum(!is.na(Y)&!is.na(X.matrix[,k]))
        model=try(summary(lm(Y~X.matrix[,k])), silent = TRUE)
        if (class(model)[1] != 'try-error') {
          b=model$coefficients[2,1]
          se=model$coefficients[2,2]
          p=model$coefficients[2,4]
          l=b-qnorm(1-sig/2)*se
          u=b+qnorm(1-sig/2)*se
          if (is.null(Adj.matrix)==F) {
            keep=!is.na(Y)&!is.na(X.matrix[,k])&apply(is.na(Adj.matrix),1,sum)==0
            n.2=sum(keep)
            y=Y[keep]
            x=X.matrix[keep,k]
            adj=Adj.matrix[keep,]
            model1=lm(y~x)
            model2=lm(y~.,data=data.frame(x,adj))        
            b1=summary(model1)$coefficients[2,1]
            se1=summary(model1)$coefficients[2,2]
            p1=summary(model1)$coefficients[2,4]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            b2=summary(model2)$coefficients[2,1]
            se2=summary(model2)$coefficients[2,2]
            p2=summary(model2)$coefficients[2,4]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
          }
        }
      }
      if (type=="logistic") {
        n.1=sum(!is.na(Y)&!is.na(X.matrix[,k]))
        model=try(summary(glm(Y~X.matrix[,k],family="binomial")), silent = TRUE)
        if (class(model)[1] != 'try-error') {
          b=model$coefficients[2,1]
          se=model$coefficients[2,2]
          p=model$coefficients[2,4]
          l=exp(b-qnorm(1-sig/2)*se)
          u=exp(b+qnorm(1-sig/2)*se)
          b=exp(b)                
          if (is.null(Adj.matrix)==F) {
            keep=!is.na(Y)&!is.na(X.matrix[,k])&apply(is.na(Adj.matrix),1,sum)==0
            n.2=sum(keep)
            y=Y[keep]
            x=X.matrix[keep,k]
            adj=Adj.matrix[keep,]
            model1=glm(y~x,family="binomial")
            model2=glm(y~.,data=data.frame(x,adj),family="binomial")        
            b1=summary(model1)$coefficients[2,1]
            se1=summary(model1)$coefficients[2,2]
            p1=summary(model1)$coefficients[2,4]
            l1=exp(b1-qnorm(1-sig/2)*se1)
            u1=exp(b1+qnorm(1-sig/2)*se1)
            b1=exp(b1)                
            b2=summary(model2)$coefficients[2,1]
            se2=summary(model2)$coefficients[2,2]
            p2=summary(model2)$coefficients[2,4]
            l2=exp(b2-qnorm(1-sig/2)*se2)
            u2=exp(b2+qnorm(1-sig/2)*se2)
            b2=exp(b2)
          }
        }
      }
      if (type=="survival") {
        n.1=sum(!is.na(Y)&!is.na(time)&!is.na(X.matrix[,k]))
        model=try(summary(coxph(Surv(time,Y)~X.matrix[,k])), silent = TRUE)
        if (class(model)[1] != 'try-error') {
          if (!is.na(model[['coefficients']][1,1])) {
            b=model$coefficients[1,1]
            se=model$coefficients[1,3]
            p=model$coefficients[1,5]
            l=exp(b-qnorm(1-sig/2)*se)
            u=exp(b+qnorm(1-sig/2)*se)
            b=exp(b)
            if (is.null(Adj.matrix)==F) {
              keep=!is.na(Y)&!is.na(time)&!is.na(X.matrix[,k])&apply(is.na(Adj.matrix),1,sum)==0
              n.2=sum(keep)
              y=Y[keep]
              t=time[keep]
              x=X.matrix[keep,k]
              adj=Adj.matrix[keep,]
              model1=coxph(Surv(t,y)~x)
              model2=coxph(Surv(t,y)~.,data=data.frame(x,adj))        
              b1=summary(model1)$coefficients[1,1]
              se1=summary(model1)$coefficients[1,3]
              p1=summary(model1)$coefficients[1,5]
              l1=exp(b1-qnorm(1-sig/2)*se1)
              u1=exp(b1+qnorm(1-sig/2)*se1)
              b1=exp(b1)                
              b2=summary(model2)$coefficients[1,1]
              se2=summary(model2)$coefficients[1,3]
              p2=summary(model2)$coefficients[1,5]
              l2=exp(b2-qnorm(1-sig/2)*se2)
              u2=exp(b2+qnorm(1-sig/2)*se2)
              b2=exp(b2)
            }
          }
        }
      }
      if (class(model)[1] != 'try-error') {
        Table[position[k],1]=paste0(.ROUND(b,digits)," (",.ROUND(l,digits),del,.ROUND(u,digits),")")
        Table[position[k],2]=.ROUND.p(p,p.digits)
        Table[position[k],3]=n.1
        if (is.null(Adj.matrix)==F) {
          Table[position[k],4]=paste0(.ROUND(b1,digits)," (",.ROUND(l1,digits),del,.ROUND(u1,digits),")")
          Table[position[k],5]=.ROUND.p(p1,p.digits)
          Table[position[k],6]=paste0(.ROUND(b2,digits)," (",.ROUND(l2,digits),del,.ROUND(u2,digits),")")
          Table[position[k],7]=.ROUND.p(p2,p.digits)
          Table[position[k],8]=n.2     
        }
      }
    } else {
      n.factor=length(levels(X.matrix[,k]))      
      if (type=="linear") {
        n.1=sum(!is.na(Y)&!is.na(X.matrix[,k]))
        model=try(summary(lm(Y~factor(X.matrix[,k]))), silent = TRUE)
        if (class(model)[1] != 'try-error') {
          b=model$coefficients[2:n.factor,1]
          se=model$coefficients[2:n.factor,2]
          p=model$coefficients[2:n.factor,4]
          l=b-qnorm(1-sig/2)*se
          u=b+qnorm(1-sig/2)*se
          overall.p=pchisq(t(b)%*%solve(vcov(model)[-1,-1])%*%b,df=n.factor-1,lower.tail=F)
          if (is.null(Adj.matrix)==F) {
            keep=!is.na(Y)&!is.na(X.matrix[,k])&apply(is.na(Adj.matrix),1,sum)==0
            n.2=sum(keep)
            y=Y[keep]
            x=X.matrix[keep,k]
            adj=Adj.matrix[keep,]
            model1=lm(y~factor(x))
            model2=lm(y~.,data=data.frame(factor(x),adj))        
            b1=summary(model1)$coefficients[2:n.factor,1]
            se1=summary(model1)$coefficients[2:n.factor,2]
            p1=summary(model1)$coefficients[2:n.factor,4]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            overall.p1=pchisq(t(b1)%*%solve(vcov(model1)[-1,-1])%*%b1,df=n.factor-1,lower.tail=F)
            b2=summary(model2)$coefficients[2:n.factor,1]
            se2=summary(model2)$coefficients[2:n.factor,2]
            p2=summary(model2)$coefficients[2:n.factor,4]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            overall.p2=pchisq(t(b2)%*%solve(vcov(model2)[2:n.factor,2:n.factor])%*%b2,df=n.factor-1,lower.tail=F)
          }
        }
      }
      if (type=="logistic") {
        n.1=sum(!is.na(Y)&!is.na(X.matrix[,k]))
        model=try(summary(glm(Y~factor(X.matrix[,k]),family="binomial")), silent = TRUE)
        if (class(model)[1] != 'try-error') {
          b=model$coefficients[2:n.factor,1]
          se=model$coefficients[2:n.factor,2]
          p=model$coefficients[2:n.factor,4]
          l=exp(b-qnorm(1-sig/2)*se)
          u=exp(b+qnorm(1-sig/2)*se)
          overall.p=pchisq(t(b)%*%solve(vcov(model)[-1,-1])%*%b,df=n.factor-1,lower.tail=F)
          b=exp(b)
          if (is.null(Adj.matrix)==F) {
            keep=!is.na(Y)&!is.na(X.matrix[,k])&apply(is.na(Adj.matrix),1,sum)==0
            n.2=sum(keep)
            y=Y[keep]
            x=X.matrix[keep,k]
            adj=Adj.matrix[keep,]
            model1=glm(y~factor(x),family="binomial")
            model2=glm(y~.,data=data.frame(factor(x),adj),family="binomial")        
            b1=summary(model1)$coefficients[2:n.factor,1]
            se1=summary(model1)$coefficients[2:n.factor,2]
            p1=summary(model1)$coefficients[2:n.factor,4]
            l1=exp(b1-qnorm(1-sig/2)*se1)
            u1=exp(b1+qnorm(1-sig/2)*se1)
            overall.p1=pchisq(t(b1)%*%solve(vcov(model1)[-1,-1])%*%b1,df=n.factor-1,lower.tail=F)
            b1=exp(b1)                
            b2=summary(model2)$coefficients[2:n.factor,1]
            se2=summary(model2)$coefficients[2:n.factor,2]
            p2=summary(model2)$coefficients[2:n.factor,4]
            l2=exp(b2-qnorm(1-sig/2)*se2)
            u2=exp(b2+qnorm(1-sig/2)*se2)
            overall.p2=pchisq(t(b2)%*%solve(vcov(model2)[2:n.factor,2:n.factor])%*%b2,df=n.factor-1,lower.tail=F)
            b2=exp(b2)
          }
        }
      }
      if (type=="survival") {
        n.1=sum(!is.na(Y)&!is.na(time)&!is.na(X.matrix[,k]))
        model=try(summary(coxph(Surv(time,Y)~factor(X.matrix[,k]))), silent = TRUE)
        if (class(model)[1] != 'try-error') {
          if (!is.na(model[['coefficients']][1,1])) {
            b=model$coefficients[1:(n.factor-1),1]
            se=model$coefficients[1:(n.factor-1),3]
            p=model$coefficients[1:(n.factor-1),5]
            l=exp(b-qnorm(1-sig/2)*se)
            u=exp(b+qnorm(1-sig/2)*se)
            overall.p=pchisq(t(b)%*%solve(vcov(coxph(Surv(time,Y)~X.matrix[,k])))%*%b,df=n.factor-1,lower.tail=F)
            b=exp(b)
            if (is.null(Adj.matrix)==F) {
              keep=!is.na(Y)&!is.na(time)&!is.na(X.matrix[,k])&apply(is.na(Adj.matrix),1,sum)==0
              n.2=sum(keep)
              y=Y[keep]
              t=time[keep]
              x=X.matrix[keep,k]
              adj=Adj.matrix[keep,]
              model1=coxph(Surv(t,y)~x)
              model2=coxph(Surv(t,y)~.,data=data.frame(factor(x),adj))        
              b1=summary(model1)$coefficients[1:(n.factor-1),1]
              se1=summary(model1)$coefficients[1:(n.factor-1),3]
              p1=summary(model1)$coefficients[1:(n.factor-1),5]
              l1=exp(b1-qnorm(1-sig/2)*se1)
              u1=exp(b1+qnorm(1-sig/2)*se1)
              overall.p1=pchisq(t(b1)%*%solve(vcov(model1))%*%b1,df=n.factor-1,lower.tail=F)
              b1=exp(b1)                
              b2=summary(model2)$coefficients[1:(n.factor-1),1]
              se2=summary(model2)$coefficients[1:(n.factor-1),3]
              p2=summary(model2)$coefficients[1:(n.factor-1),5]
              l2=exp(b2-qnorm(1-sig/2)*se2)
              u2=exp(b2+qnorm(1-sig/2)*se2)
              overall.p2=pchisq(t(b2)%*%solve(vcov(model2)[1:(n.factor-1),1:(n.factor-1)])%*%b2,df=n.factor-1,lower.tail=F)
              b2=exp(b2)
            }
          }
        }
      }
      if (class(model)[1] != 'try-error') {
        if (!is.na(model[['coefficients']][1,1])) {
          Table[(position[k]-n.factor+2):position[k],1]=paste0(.ROUND(b,digits)," (",.ROUND(l,digits),del,.ROUND(u,digits),")")
          Table[(position[k]-n.factor+2):position[k],2]=.ROUND.p(p,p.digits)
          if (type=="linear") {Table[(position[k]-n.factor+1),1]=.ROUND(0,digits)} else {
            Table[(position[k]-n.factor+1),1]=.ROUND(1,digits)
          }
          Table[(position[k]-n.factor),2]=.ROUND.p(overall.p,p.digits)
          Table[(position[k]-n.factor),3]=n.1      
          if (is.null(Adj.matrix)==F) {
            Table[(position[k]-n.factor+2):position[k],4]=paste0(.ROUND(b1,digits)," (",.ROUND(l1,digits),del,.ROUND(u1,digits),")")
            Table[(position[k]-n.factor+2):position[k],5]=.ROUND.p(p1,p.digits)
            Table[(position[k]-n.factor+2):position[k],6]=paste0(.ROUND(b2,digits)," (",.ROUND(l2,digits),del,.ROUND(u2,digits),")")
            Table[(position[k]-n.factor+2):position[k],7]=.ROUND.p(p2,p.digits) 
            if (type=="linear") {
              Table[(position[k]-n.factor+1),4]=.ROUND(0,digits)
              Table[(position[k]-n.factor+1),6]=.ROUND(0,digits)
            } else {
              Table[(position[k]-n.factor+1),4]=.ROUND(1,digits)
              Table[(position[k]-n.factor+1),6]=.ROUND(1,digits)
            }
            Table[(position[k]-n.factor),5]=.ROUND.p(overall.p1,p.digits)            
            Table[(position[k]-n.factor),7]=.ROUND.p(overall.p2,p.digits)
            Table[(position[k]-n.factor),8]=n.2     
          }
        }
      }
    }
  }
  Table.Final=cbind(rownames(Table),Table)
  colnames(Table.Final)[1]="Independent variable"
  rownames(Table.Final)=1:nrow(Table)
  return(Table.Final)
}

Table3 <- function(X,Y,M.matrix,Adj.matrix=NULL,time=NULL,type,sig=0.05,digits=2,p.digits=3,Factor=NULL,del=" to ") {
  if ("Table"=="Table") {
    #if (length(M.matrix)==length(Y)) {M.matrix=matrix(M.matrix,length(Y),1)}
    if ((class(M.matrix)!="matrix"|class(M.matrix)!="data.frame")&length(M.matrix)==length(Y)){M.matrix=matrix(M.matrix,length(Y),1)}
    if (length(Adj.matrix)==length(Y)) {Adj.matrix=matrix(Adj.matrix,length(Y),1)}  
    if (is.null(Factor)) {if (is.factor(X[,1])) {Factor=1} else {Factor=0}}
    if (Factor==1) {if (table(X[,1])[1]<table(X[,1])[length(levels(X[,1]))]) {X[,1]=factor(X[,1], levels = levels(X[,1])[length(levels(X[,1])):1])}}
    
    n.row=numeric()
    row.names=NULL
    row.names2=NULL
    for (k in 1:ncol(M.matrix)) {
      if (Factor==0) {
        n.row=c(n.row,1+length(levels(M.matrix[,k])))
        row.names=c(row.names,colnames(M.matrix)[k],paste0(colnames(M.matrix)[k],":",levels(factor(M.matrix[,k]))))
      } else {
        n.row=c(n.row,1+length(levels(M.matrix[,k]))*length(levels(X[,1])))
        row.names=c(row.names,colnames(M.matrix)[k])
        row.names2=c(row.names2,"")
        for (j in 1:length(levels(M.matrix[,k]))) {
          row.names=c(row.names,paste0(colnames(M.matrix)[k],":",c(levels(factor(M.matrix[,k]))[j])))
          row.names=c(row.names,rep("",length(levels(X[,1]))-1))
          row.names2=c(row.names2,paste0(names(X),":",levels(factor(X[,1]))))
        }
      }
    }
    
    if (is.null(Adj.matrix)) {
      if (is.null(row.names2)) {
        Table=matrix("",nrow=sum(n.row),ncol=4)
        if (type=="linear") {coef.name="beta"}
        if (type=="logistic") {coef.name="OR"}
        if (type=="survival") {coef.name="HR"}  
        colnames(Table)=c(paste0(coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)","n")
        rownames(Table)=row.names
        position=cumsum(n.row)
      } else {
        Table=matrix("",nrow=sum(n.row),ncol=5)
        if (type=="linear") {coef.name="beta"}
        if (type=="logistic") {coef.name="OR"}
        if (type=="survival") {coef.name="HR"}  
        colnames(Table)=c(names(X),paste0(coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)","n")
        rownames(Table)=row.names
        Table[,1]=row.names2
        position=cumsum(n.row)
      }
    } else {
      if (is.null(row.names2)) {
        Table=matrix("",nrow=sum(n.row),ncol=11)
        if (type=="linear") {coef.name="beta"}
        if (type=="logistic") {coef.name="OR"}
        if (type=="survival") {coef.name="HR"}  
        colnames(Table)=c(paste0("Subset 1 Crude-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)","n",paste0("Subset 2 Crude-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)",paste0("Subset 2 Adj-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)","n")
        rownames(Table)=row.names
        position=cumsum(n.row)
      } else {
        Table=matrix("",nrow=sum(n.row),ncol=12)
        if (type=="linear") {coef.name="beta"}
        if (type=="logistic") {coef.name="OR"}
        if (type=="survival") {coef.name="HR"}  
        colnames(Table)=c(names(X),paste0("Subset 1 Crude-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)","n",paste0("Subset 2 Crude-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)",paste0("Subset 2 Adj-",coef.name," (",round((1-sig)*100,2),"% CI)"),"p-value","p-value(interaction)","n")
        rownames(Table)=row.names
        Table[,1]=row.names2
        position=cumsum(n.row)
      }
    }
  }
  
  #######################################
  
  for (k in 1:ncol(M.matrix)) {
    if (type=="linear") {
      if (1==1) {
        n.1=tapply(!is.na(Y)&!is.na(X[,1]),M.matrix[,k],sum)
        for (j in 1:length(levels(M.matrix[,k]))) {
          tn=table(X[!is.na(Y)&!is.na(X[,1])&M.matrix[,k]==levels(M.matrix[,k])[j],1])
          if (sum(tn==0)==0) {
            model=summary(lm(Y[M.matrix[,k]==levels(M.matrix[,k])[j]]~X[M.matrix[,k]==levels(M.matrix[,k])[j],1]))
            b=model$coefficients[-1,1]
            se=model$coefficients[-1,2]
            p=model$coefficients[-1,4]
            l=b-qnorm(1-sig/2)*se
            u=b+qnorm(1-sig/2)*se
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j),1]=paste0(.ROUND(b,digits)," (",.ROUND(l,digits),del,.ROUND(u,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j),2]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j),4]=n.1[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,2]=.ROUND(0,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),2]=paste0(.ROUND(b,digits)," (",.ROUND(l,digits),del,.ROUND(u,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),3]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,5]=n.1[j]
            }
          } else if (sum(tn!=0)>=2&tn[1]>0) {
            model=summary(lm(Y[M.matrix[,k]==levels(M.matrix[,k])[j]]~X[M.matrix[,k]==levels(M.matrix[,k])[j],1]))
            b=rep(NA,rep(length(levels(X[,1])))-1);b[(tn!=0)[-1]]=model$coefficients[-1,1]
            se=rep(NA,rep(length(levels(X[,1])))-1);se[(tn!=0)[-1]]=model$coefficients[-1,2]
            p=rep(NA,rep(length(levels(X[,1])))-1);p[(tn!=0)[-1]]=model$coefficients[-1,4]
            l=b-qnorm(1-sig/2)*se
            u=b+qnorm(1-sig/2)*se
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j),1]=paste0(.ROUND(b,digits)," (",.ROUND(l,digits),del,.ROUND(u,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j),2]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j),4]=n.1[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,2]=.ROUND(0,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),2]=paste0(.ROUND(b,digits)," (",.ROUND(l,digits),del,.ROUND(u,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),3]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,5]=n.1[j]
            }
          }
        }
        model=summary(lm(Y~X[,1]+M.matrix[,k]+X[,1]*M.matrix[,k]))
        if (is.null(row.names2)) {
          b=model$coefficients[-c(1:(1+length(levels(M.matrix[,k])))),1]
          v=vcov(model)[-c(1:(1+length(levels(M.matrix[,k])))),-c(1:(1+length(levels(M.matrix[,k]))))]
          v=as.matrix(v)
          v=v[!is.na(b),!is.na(b)]
          b=b[!is.na(b)]
          overall.p=pchisq(t(b)%*%solve(v)%*%b,df=length(b),lower.tail=F)
          Table[position[k]-n.row[k]+1,3]=.ROUND.p(overall.p,p.digits)
        } else {
          b=model$coefficients[-c(1:(length(levels(X[,1]))+length(levels(M.matrix[,k]))-1)),1]
          v=vcov(model)[-c(1:(length(levels(X[,1]))+length(levels(M.matrix[,k]))-1)),-c(1:(length(levels(X[,1]))+length(levels(M.matrix[,k]))-1))]
          v=as.matrix(v)
          v=v[!is.na(b),!is.na(b)]
          b=b[!is.na(b)]
          overall.p=pchisq(t(b)%*%solve(v)%*%b,df=length(b),lower.tail=F)
          Table[position[k]-n.row[k]+1,4]=.ROUND.p(overall.p,p.digits)
        }
      }
      if (!is.null(Adj.matrix)) {
        keep=!is.na(Y)&!is.na(X[,1])&apply(is.na(Adj.matrix),1,sum)==0
        y=Y[keep]
        x=X[keep,1];x=data.frame(x)
        if (Factor==1) {x[,1]=factor(x[,1],levels = levels(X[,1]))}
        adj=Adj.matrix[keep,]
        m=M.matrix[keep,k]
        n.2=tapply(!is.na(y)&!is.na(x[,1]),m,sum)
        for (j in 1:length(levels(M.matrix[,k]))) {
          tn=table(x[!is.na(y)&!is.na(x[,1])&m==levels(M.matrix[,k])[j],1])
          if (sum(tn==0)==0) {
            model1=summary(lm(y[m==levels(m)[j]]~x[m==levels(m)[j],1]))
            new.data=data.frame(x,adj)[m==levels(m)[j],]
            save.variable=apply(new.data,2,function(x) {return(sd(x,na.rm=T))})!=0
            save.variable[is.na(save.variable)]=apply(data.frame(new.data[,is.na(save.variable)]),2,function(x) {length(levels(factor(x)))})>1
            new.data=new.data[,save.variable]
            model2=summary(lm(y[m==levels(m)[j]]~.,data=data.frame(new.data)))
            b1=model1$coefficients[-1,1]
            se1=model1$coefficients[-1,2]
            p1=model1$coefficients[-1,4]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            lvls=length(levels(x[,1]))
            if (lvls==0) {lvls=2}
            b2=model2$coefficients[2:lvls,1]
            se2=model2$coefficients[2:lvls,2]
            p2=model2$coefficients[2:lvls,4]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,5]=paste0(.ROUND(b1,digits)," (",.ROUND(l1,digits),del,.ROUND(u1,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,6]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,8]=paste0(.ROUND(b2,digits)," (",.ROUND(l2,digits),del,.ROUND(u2,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,9]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,11]=n.2[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,6]=.ROUND(0,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),6]=paste0(.ROUND(b1,digits)," (",.ROUND(l1,digits),del,.ROUND(u1,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),7]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,9]=.ROUND(0,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),9]=paste0(.ROUND(b2,digits)," (",.ROUND(l2,digits),del,.ROUND(u2,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),10]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,12]=n.2[j]
            }
          } else if (sum(tn!=0)>=2&tn[1]>0) {
            model1=summary(lm(y[m==levels(m)[j]]~x[m==levels(m)[j],1]))
            new.data=data.frame(x,adj)[m==levels(m)[j],]
            save.variable=apply(new.data,2,function(x) {return(sd(x,na.rm=T))})!=0
            save.variable[is.na(save.variable)]=apply(data.frame(new.data[,is.na(save.variable)]),2,function(x) {length(levels(factor(x)))})>1
            new.data=new.data[,save.variable]
            model2=summary(lm(y[m==levels(m)[j]]~.,data=data.frame(new.data)))
            b1=rep(NA,rep(length(levels(x[,1])))-1);b1[(tn!=0)[-1]]=model1$coefficients[-1,1]
            se1=rep(NA,rep(length(levels(x[,1])))-1);se1[(tn!=0)[-1]]=model1$coefficients[-1,2]
            p1=rep(NA,rep(length(levels(x[,1])))-1);p1[(tn!=0)[-1]]=model1$coefficients[-1,4]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            lvls=length(levels(x[,1]))
            if (lvls==0) {lvls=2}
            b2=rep(NA,rep(length(levels(x[,1])))-1);b2[(tn!=0)[-1]]=model2$coefficients[2:(lvls-sum((tn!=0)[-1])),1]
            se2=rep(NA,rep(length(levels(x[,1])))-1);se2[(tn!=0)[-1]]=model2$coefficients[2:(lvls-sum((tn!=0)[-1])),2]
            p2=rep(NA,rep(length(levels(x[,1])))-1);p2[(tn!=0)[-1]]=model2$coefficients[2:(lvls-sum((tn!=0)[-1])),4]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,5]=paste0(.ROUND(b1,digits)," (",.ROUND(l1,digits),del,.ROUND(u1,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,6]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,8]=paste0(.ROUND(b2,digits)," (",.ROUND(l2,digits),del,.ROUND(u2,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,9]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,11]=n.2[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,6]=.ROUND(0,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),6]=paste0(.ROUND(b1,digits)," (",.ROUND(l1,digits),del,.ROUND(u1,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),7]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,9]=.ROUND(0,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),9]=paste0(.ROUND(b2,digits)," (",.ROUND(l2,digits),del,.ROUND(u2,digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),10]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,12]=n.2[j]
            }
          }
        }
        model1=summary(lm(y~x[,1]+m+x[,1]*m))
        model2=summary(lm(y~x[,1]+m+x[,1]*m+.,data=data.frame(adj)))
        if (is.null(row.names2)) {
          b1=model1$coefficients[(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients),1]
          v1=vcov(model1)[(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients),(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients)]
          v1=as.matrix(v1)
          v1=v1[!is.na(b1),!is.na(b1)]
          b1=b1[!is.na(b1)]
          overall.p1=pchisq(t(b1)%*%solve(v1)%*%b1,df=length(b1),lower.tail=F)
          b2=model2$coefficients[(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients),1]
          v2=vcov(model2)[(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients),(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients)]
          v2=as.matrix(v2)
          v2=v2[!is.na(b2),!is.na(b2)]
          b2=b2[!is.na(b2)]
          overall.p2=pchisq(t(b2)%*%solve(v2)%*%b2,df=length(b2),lower.tail=F)
          Table[position[k]-n.row[k]+1,7]=.ROUND.p(overall.p1,p.digits)
          Table[position[k]-n.row[k]+1,10]=.ROUND.p(overall.p2,p.digits)
        } else {
          b1=model1$coefficients[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)-1)),1]
          v1=vcov(model1)[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)-1)),-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)-1))]
          v1=as.matrix(v1)
          v1=v1[!is.na(b1),!is.na(b1)]
          b1=b1[!is.na(b1)]
          overall.p1=pchisq(t(b1)%*%solve(v1)%*%b1,df=length(b1),lower.tail=F)
          n.df=NULL
          for (i in 1:ncol(data.frame(adj))) {
            n.df[i]=length(levels(data.frame(adj)[,i]))-1
            if (n.df[i]==-1) {n.df[i]=1}
            if (all.equal(factor(data.frame(adj)[,i]),factor(m))[1]==TRUE) {n.df[i]=0}
            if (all.equal(factor(data.frame(adj)[,i]),factor(x[,1]))[1]==TRUE) {n.df[i]=0}
          }
          b2=model2$coefficients[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)+sum(n.df)-1)),1]
          v2=vcov(model2)[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)+sum(n.df)-1)),-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)+sum(n.df)-1))]
          v2=as.matrix(v2)
          v2=v2[!is.na(b2),!is.na(b2)]
          b2=b2[!is.na(b2)]
          overall.p2=pchisq(t(b2)%*%solve(v2)%*%b2,df=length(b2),lower.tail=F)
          Table[position[k]-n.row[k]+1,8]=.ROUND.p(overall.p1,p.digits)
          Table[position[k]-n.row[k]+1,11]=.ROUND.p(overall.p2,p.digits)
        }
      }
    }
    if (type=="logistic") {
      if (1==1) {
        n.1=tapply(!is.na(Y)&!is.na(X[,1]),M.matrix[,k],sum)
        for (j in 1:length(levels(M.matrix[,k]))) {
          tn=table(X[!is.na(Y)&!is.na(X[,1])&M.matrix[,k]==levels(M.matrix[,k])[j],1])
          if (sum(tn==0)==0) {
            model=summary(glm(Y[M.matrix[,k]==levels(M.matrix[,k])[j]]~X[M.matrix[,k]==levels(M.matrix[,k])[j],1],family="binomial"))
            b=model$coefficients[-1,1]
            se=model$coefficients[-1,2]
            p=model$coefficients[-1,4]
            l=b-qnorm(1-sig/2)*se
            u=b+qnorm(1-sig/2)*se
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j),1]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j),2]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j),4]=n.1[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,2]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),2]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),3]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,5]=n.1[j]
            }
          } else if (sum(tn!=0)>=2&tn[1]>0) {
            model=summary(glm(Y[M.matrix[,k]==levels(M.matrix[,k])[j]]~X[M.matrix[,k]==levels(M.matrix[,k])[j],1],family="binomial"))
            b=rep(NA,rep(length(levels(X[,1])))-1);b[(tn!=0)[-1]]=model$coefficients[-1,1]
            se=rep(NA,rep(length(levels(X[,1])))-1);se[(tn!=0)[-1]]=model$coefficients[-1,2]
            p=rep(NA,rep(length(levels(X[,1])))-1);p[(tn!=0)[-1]]=model$coefficients[-1,4]
            l=b-qnorm(1-sig/2)*se
            u=b+qnorm(1-sig/2)*se
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j),1]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j),2]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j),4]=n.1[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,2]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),2]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),3]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,5]=n.1[j]
            }         
          }
        }
        model=summary(glm(Y~X[,1]+M.matrix[,k]+X[,1]*M.matrix[,k],family="binomial"))
        if (is.null(row.names2)) {
          b=model$coefficients[-c(1:(1+length(levels(M.matrix[,k])))),1]
          v=vcov(model)[-c(1:(1+length(levels(M.matrix[,k])))),-c(1:(1+length(levels(M.matrix[,k]))))]
          v=as.matrix(v)
          v=v[!is.na(b),!is.na(b)]
          b=b[!is.na(b)]
          overall.p=pchisq(t(b)%*%solve(v)%*%b,df=length(b),lower.tail=F)
          Table[position[k]-n.row[k]+1,3]=.ROUND.p(overall.p,p.digits)
        } else {
          b=model$coefficients[-c(1:(length(levels(X[,1]))+length(levels(M.matrix[,k]))-1)),1]
          v=vcov(model)[-c(1:(length(levels(X[,1]))+length(levels(M.matrix[,k]))-1)),-c(1:(length(levels(X[,1]))+length(levels(M.matrix[,k]))-1))]
          v=as.matrix(v)
          v=v[!is.na(b),!is.na(b)]
          b=b[!is.na(b)]
          overall.p=pchisq(t(b)%*%solve(v)%*%b,df=length(b),lower.tail=F)
          Table[position[k]-n.row[k]+1,4]=.ROUND.p(overall.p,p.digits)
        }
      }
      if (!is.null(Adj.matrix)) {
        keep=!is.na(Y)&!is.na(X[,1])&apply(is.na(Adj.matrix),1,sum)==0
        y=Y[keep]
        x=X[keep,1];x=data.frame(x)
        if (Factor==1) {x[,1]=factor(x[,1],levels = levels(X[,1]))}
        adj=Adj.matrix[keep,]
        m=M.matrix[keep,k]
        n.2=tapply(!is.na(y)&!is.na(x[,1]),m,sum)
        for (j in 1:length(levels(M.matrix[,k]))) {
          tn=table(x[!is.na(y)&!is.na(x[,1])&m==levels(M.matrix[,k])[j],1])
          if (sum(tn==0)==0) {
            model1=summary(glm(y[m==levels(m)[j]]~x[m==levels(m)[j],1],family="binomial"))
            new.data=data.frame(x,adj)[m==levels(m)[j],]
            save.variable=apply(new.data,2,function(x) {return(sd(x,na.rm=T))})!=0
            save.variable[is.na(save.variable)]=apply(data.frame(new.data[,is.na(save.variable)]),2,function(x) {length(levels(factor(x)))})>1
            new.data=new.data[,save.variable]
            model2=summary(glm(y[m==levels(m)[j]]~.,data=data.frame(new.data),family="binomial"))
            b1=model1$coefficients[-1,1]
            se1=model1$coefficients[-1,2]
            p1=model1$coefficients[-1,4]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            lvls=length(levels(x[,1]))
            if (lvls==0) {lvls=2}
            b2=model2$coefficients[2:lvls,1]
            se2=model2$coefficients[2:lvls,2]
            p2=model2$coefficients[2:lvls,4]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,5]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,6]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,8]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,9]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,11]=n.2[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,6]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),6]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),7]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,9]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),9]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),10]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,12]=n.2[j]
            }
          } else if (sum(tn!=0)>=2&tn[1]>0) {
            model1=summary(glm(y[m==levels(m)[j]]~x[m==levels(m)[j],1],family="binomial"))
            new.data=data.frame(x,adj)[m==levels(m)[j],]
            save.variable=apply(new.data,2,function(x) {return(sd(x,na.rm=T))})!=0
            save.variable[is.na(save.variable)]=apply(data.frame(new.data[,is.na(save.variable)]),2,function(x) {length(levels(factor(x)))})>1
            new.data=new.data[,save.variable]
            model2=summary(glm(y[m==levels(m)[j]]~.,data=data.frame(new.data),family="binomial"))
            b1=rep(NA,rep(length(levels(x[,1])))-1);b1[(tn!=0)[-1]]=model1$coefficients[-1,1]
            se1=rep(NA,rep(length(levels(x[,1])))-1);se1[(tn!=0)[-1]]=model1$coefficients[-1,2]
            p1=rep(NA,rep(length(levels(x[,1])))-1);p1[(tn!=0)[-1]]=model1$coefficients[-1,4]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            lvls=length(levels(x[,1]))
            if (lvls==0) {lvls=2}
            b2=rep(NA,rep(length(levels(x[,1])))-1);b2[(tn!=0)[-1]]=model2$coefficients[2:(lvls-sum((tn!=0)[-1])),1]
            se2=rep(NA,rep(length(levels(x[,1])))-1);se2[(tn!=0)[-1]]=model2$coefficients[2:(lvls-sum((tn!=0)[-1])),2]
            p2=rep(NA,rep(length(levels(x[,1])))-1);p2[(tn!=0)[-1]]=model2$coefficients[2:(lvls-sum((tn!=0)[-1])),4]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,5]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,6]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,8]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,9]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,11]=n.2[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,6]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),6]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),7]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,9]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),9]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),10]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,12]=n.2[j]
            }  
          }
        }
        model1=summary(glm(y~x[,1]+m+x[,1]*m,family="binomial"))
        model2=summary(glm(y~x[,1]+m+x[,1]*m+.,data=data.frame(adj),family="binomial"))
        if (is.null(row.names2)) {
          b1=model1$coefficients[(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients),1]
          v1=vcov(model1)[(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients),(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients)]
          v1=as.matrix(v1)
          v1=v1[!is.na(b1),!is.na(b1)]
          b1=b1[!is.na(b1)]
          overall.p1=pchisq(t(b1)%*%solve(v1)%*%b1,df=length(b1),lower.tail=F)
          b2=model2$coefficients[(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients),1]
          v2=vcov(model2)[(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients),(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients)]
          v2=as.matrix(v2)
          v2=v2[!is.na(b2),!is.na(b2)]
          b2=b2[!is.na(b2)]
          overall.p2=pchisq(t(b2)%*%solve(v2)%*%b2,df=length(b2),lower.tail=F)
          Table[position[k]-n.row[k]+1,7]=.ROUND.p(overall.p1,p.digits)
          Table[position[k]-n.row[k]+1,10]=.ROUND.p(overall.p2,p.digits)
        } else {
          b1=model1$coefficients[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)-1)),1]
          v1=vcov(model1)[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)-1)),-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)-1))]
          v1=as.matrix(v1)
          v1=v1[!is.na(b1),!is.na(b1)]
          b1=b1[!is.na(b1)]
          overall.p1=pchisq(t(b1)%*%solve(v1)%*%b1,df=length(b1),lower.tail=F)
          n.df=NULL
          for (i in 1:ncol(data.frame(adj))) {
            n.df[i]=length(levels(data.frame(adj)[,i]))-1
            if (n.df[i]==-1) {n.df[i]=1}
            if (all.equal(factor(data.frame(adj)[,i]),factor(m))[1]==TRUE) {n.df[i]=0}
            if (all.equal(factor(data.frame(adj)[,i]),factor(x[,1]))[1]==TRUE) {n.df[i]=0}
          }
          b2=model2$coefficients[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)+sum(n.df)-1)),1]
          v2=vcov(model2)[-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)+sum(n.df)-1)),-c(1:(sum(table(x[,1])!=0)+sum(table(m)!=0)+sum(n.df)-1))]
          v2=as.matrix(v2)
          v2=v2[!is.na(b2),!is.na(b2)]
          b2=b2[!is.na(b2)]
          overall.p2=pchisq(t(b2)%*%solve(v2)%*%b2,df=length(b2),lower.tail=F)
          Table[position[k]-n.row[k]+1,8]=.ROUND.p(overall.p1,p.digits)
          Table[position[k]-n.row[k]+1,11]=.ROUND.p(overall.p2,p.digits)
        }
      }
    }
    if (type=="survival") {
      if (1==1) {
        n.1=tapply(!is.na(Y)&!is.na(time)&!is.na(X[,1]),M.matrix[,k],sum)
        for (j in 1:length(levels(M.matrix[,k]))) {
          tn=table(X[!is.na(Y)&!is.na(time)&!is.na(X[,1])&M.matrix[,k]==levels(M.matrix[,k])[j],1])
          if (sum(tn==0)==0) {
            model=summary(coxph(Surv(time[M.matrix[,k]==levels(M.matrix[,k])[j]],Y[M.matrix[,k]==levels(M.matrix[,k])[j]])~X[M.matrix[,k]==levels(M.matrix[,k])[j],1]))
            b=model$coefficients[,1]
            se=model$coefficients[,3]
            p=model$coefficients[,5]
            l=b-qnorm(1-sig/2)*se
            u=b+qnorm(1-sig/2)*se
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j),1]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j),2]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j),4]=n.1[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,2]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),2]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),3]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,5]=n.1[j]
            }
          } else if (sum(tn!=0)>=2&tn[1]>0) {
            model=summary(coxph(Surv(time[M.matrix[,k]==levels(M.matrix[,k])[j]],Y[M.matrix[,k]==levels(M.matrix[,k])[j]])~X[M.matrix[,k]==levels(M.matrix[,k])[j],1]))
            b=model$coefficients[,1]
            se=model$coefficients[,3]
            p=model$coefficients[,5]
            l=b-qnorm(1-sig/2)*se
            u=b+qnorm(1-sig/2)*se
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j),1]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j),2]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j),4]=n.1[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,2]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),2]=paste0(.ROUND(exp(b),digits)," (",.ROUND(exp(l),digits),del,.ROUND(exp(u),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),3]=.ROUND.p(p,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,5]=n.1[j]
            }
          }
        }
        model=summary(coxph(Surv(time,Y)~X[,1]+M.matrix[,k]+X[,1]*M.matrix[,k]))
        if (is.null(row.names2)) {
          b=model$coefficients[(nrow(model$coefficients)-(length(levels(M.matrix[,k]))-1)+1):nrow(model$coefficients),1]
          v=vcov(coxph(Surv(time,Y)~X[,1]+M.matrix[,k]+X[,1]*M.matrix[,k]))[(nrow(model$coefficients)-(length(levels(M.matrix[,k]))-1)+1):nrow(model$coefficients),(nrow(model$coefficients)-(length(levels(M.matrix[,k]))-1)+1):nrow(model$coefficients)]
          v=as.matrix(v)
          v=v[!is.na(b),!is.na(b)]
          b=b[!is.na(b)]
          overall.p=pchisq(t(b)%*%solve(v)%*%b,df=length(b),lower.tail=F)
          Table[position[k]-n.row[k]+1,3]=.ROUND.p(overall.p,p.digits)
        } else {
          b=model$coefficients[(nrow(model$coefficients)-(length(levels(X[,1]))-1)*(length(levels(M.matrix[,k]))-1)+1):nrow(model$coefficients),1]
          v=vcov(coxph(Surv(time,Y)~X[,1]+M.matrix[,k]+X[,1]*M.matrix[,k]))[(nrow(model$coefficients)-(length(levels(X[,1]))-1)*(length(levels(M.matrix[,k]))-1)+1):nrow(model$coefficients),(nrow(model$coefficients)-(length(levels(X[,1]))-1)*(length(levels(M.matrix[,k]))-1)+1):nrow(model$coefficients)]
          v=as.matrix(v)
          v=v[!is.na(b),!is.na(b)]
          b=b[!is.na(b)]
          overall.p=pchisq(t(b)%*%solve(v)%*%b,df=length(b),lower.tail=F)
          Table[position[k]-n.row[k]+1,4]=.ROUND.p(overall.p,p.digits)
        }
      }
      if (!is.null(Adj.matrix)) {
        keep=!is.na(Y)&!is.na(time)&!is.na(X[,1])&apply(is.na(Adj.matrix),1,sum)==0
        y=Y[keep]
        t=time[keep]
        x=X[keep,1];x=data.frame(x)
        if (Factor==1) {x[,1]=factor(x[,1],levels = levels(X[,1]))}
        adj=Adj.matrix[keep,]
        m=M.matrix[keep,k]
        n.2=tapply(!is.na(y)&!is.na(t)&!is.na(x[,1]),m,sum)
        for (j in 1:length(levels(M.matrix[,k]))) {
          tn=table(x[!is.na(y)&!is.na(x[,1])&m==levels(M.matrix[,k])[j],1])
          if (sum(tn==0)==0) {
            model1=summary(coxph(Surv(t[m==levels(m)[j]],y[m==levels(m)[j]])~x[m==levels(m)[j],1]))
            new.data=data.frame(x,adj)[m==levels(m)[j],]
            save.variable=apply(new.data,2,function(x) {return(sd(x,na.rm=T))})!=0
            save.variable[is.na(save.variable)]=apply(data.frame(new.data[,is.na(save.variable)]),2,function(x) {length(levels(factor(x)))})>1
            new.data=new.data[,save.variable]
            model2=summary(coxph(Surv(t[m==levels(m)[j]],y[m==levels(m)[j]])~.,data=data.frame(new.data)))
            b1=model1$coefficients[,1]
            se1=model1$coefficients[,3]
            p1=model1$coefficients[,5]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            lvls=length(levels(x[,1]))
            if (lvls==0) {lvls=2}
            b2=model2$coefficients[1:(lvls-1),1]
            se2=model2$coefficients[1:(lvls-1),3]
            p2=model2$coefficients[1:(lvls-1),5]
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,5]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,6]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,8]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,9]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,11]=n.2[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,6]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),6]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),7]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,9]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),9]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),10]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,12]=n.2[j]
            }
          } else if (sum(tn!=0)>=2&tn[1]>0) {
            model1=summary(coxph(Surv(t[m==levels(m)[j]],y[m==levels(m)[j]])~x[m==levels(m)[j],1]))
            new.data=data.frame(x,adj)[m==levels(m)[j],]
            save.variable=apply(new.data,2,function(x) {return(sd(x,na.rm=T))})!=0
            save.variable[is.na(save.variable)]=apply(data.frame(new.data[,is.na(save.variable)]),2,function(x) {length(levels(factor(x)))})>1
            new.data=new.data[,save.variable]
            model2=summary(coxph(Surv(t[m==levels(m)[j]],y[m==levels(m)[j]])~.,data=data.frame(new.data)))
            b1=model1$coefficients[,1]
            se1=model1$coefficients[,3]
            p1=model1$coefficients[,5]
            l1=b1-qnorm(1-sig/2)*se1
            u1=b1+qnorm(1-sig/2)*se1
            lvls=length(levels(x[,1]))
            if (lvls==0) {lvls=2}
            b2=model2$coefficients[1:(lvls-1),1]
            se2=model2$coefficients[1:(lvls-1),3]
            p2=model2$coefficients[1:(lvls-1),5]
            b2[is.na(p2)]<-NA
            se2[is.na(p2)]<-NA
            l2=b2-qnorm(1-sig/2)*se2
            u2=b2+qnorm(1-sig/2)*se2
            if (is.null(row.names2)) {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,5]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,6]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,8]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,9]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)+1,11]=n.2[j]
            } else {
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,6]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),6]=paste0(.ROUND(exp(b1),digits)," (",.ROUND(exp(l1),digits),del,.ROUND(exp(u1),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),7]=.ROUND.p(p1,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,9]=.ROUND(1,digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),9]=paste0(.ROUND(exp(b2),digits)," (",.ROUND(exp(l2),digits),del,.ROUND(exp(u2),digits),")")
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+(2:length(levels(X[,1]))),10]=.ROUND.p(p2,p.digits)
              Table[position[k]-(length(levels(M.matrix[,k]))-j+1)*length(levels(X[,1]))+1,12]=n.2[j]
            }
          }
        }
        model1=summary(coxph(Surv(t,y)~x[,1]+m+x[,1]*m))
        model2=summary(coxph(Surv(t,y)~x[,1]+m+x[,1]*m+.,data=data.frame(adj)))
        if (is.null(row.names2)) {
          b1=model1$coefficients[(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients),1]
          v1=vcov(coxph(Surv(t,y)~x[,1]+m+x[,1]*m))[(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients),(nrow(model1$coefficients)-(length(levels(m))-1)+1):nrow(model1$coefficients)]
          v1=as.matrix(v1)
          v1=v1[!is.na(b1),!is.na(b1)]
          b1=b1[!is.na(b1)]
          overall.p1=pchisq(t(b1)%*%solve(v1)%*%b1,df=length(b1),lower.tail=F)
          b2=model2$coefficients[(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients),1]
          v2=vcov(coxph(Surv(t,y)~x[,1]+m+x[,1]*m+.,data=data.frame(adj)))[(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients),(nrow(model2$coefficients)-(length(levels(m))-1)+1):nrow(model2$coefficients)]
          v2=as.matrix(v2)
          v2=v2[!is.na(b2),!is.na(b2)]
          b2=b2[!is.na(b2)]
          overall.p2=pchisq(t(b2)%*%solve(v2)%*%b2,df=length(b2),lower.tail=F)
          Table[position[k]-n.row[k]+1,7]=.ROUND.p(overall.p1,p.digits)
          Table[position[k]-n.row[k]+1,10]=.ROUND.p(overall.p2,p.digits)
        } else {
          b1=model1$coefficients[(nrow(model1$coefficients)-(length(levels(x[,1]))-1)*(length(levels(m))-1)+1):nrow(model1$coefficients),1]
          v1=vcov(coxph(Surv(t,y)~x[,1]+m+x[,1]*m))[(nrow(model1$coefficients)-(length(levels(x[,1]))-1)*(length(levels(m))-1)+1):nrow(model1$coefficients),(nrow(model1$coefficients)-(length(levels(x[,1]))-1)*(length(levels(m))-1)+1):nrow(model1$coefficients)]
          v1=as.matrix(v1)
          v1=v1[!is.na(b1),!is.na(b1)]
          b1=b1[!is.na(b1)]
          overall.p1=pchisq(t(b1)%*%solve(v1)%*%b1,df=length(b1),lower.tail=F)
          b2=model2$coefficients[(nrow(model2$coefficients)-(length(levels(x[,1]))-1)*(length(levels(m))-1)+1):nrow(model2$coefficients),1]
          v2=vcov(coxph(Surv(t,y)~x[,1]+m+x[,1]*m+.,data=data.frame(adj)))[(nrow(model2$coefficients)-(length(levels(x[,1]))-1)*(length(levels(m))-1)+1):nrow(model2$coefficients),(nrow(model2$coefficients)-(length(levels(x[,1]))-1)*(length(levels(m))-1)+1):nrow(model2$coefficients)]
          v2=as.matrix(v2)
          v2=v2[!is.na(b2),!is.na(b2)]
          b2=b2[!is.na(b2)]
          overall.p2=pchisq(t(b2)%*%solve(v2)%*%b2,df=length(b2),lower.tail=F)
          Table[position[k]-n.row[k]+1,8]=.ROUND.p(overall.p1,p.digits)
          Table[position[k]-n.row[k]+1,11]=.ROUND.p(overall.p2,p.digits)
        }
      }
    }
  }
  Table.Final=cbind(rownames(Table),Table)
  colnames(Table.Final)[1]="Stratified variable"
  rownames(Table.Final)=1:nrow(Table)
  return(Table.Final)
}

Table2doc <- function (table_list, filename = 'test.doc') {
  
  rtffile <- RTF(filename, width = 8.5, height = 11,font.size = 10, omi = c(0.5, 0.5, 0.5, 0.5))
  
  for (j in 1:length(table_list)) {
    
    addHeader(rtffile, paste0("Table ", j))
    ALIGN <- c("L", rep("C", ncol(table_list[[j]]) - 1))
    new_table <- table_list[[j]]
    addTable(rtffile, new_table, NA.string = "", row.names = FALSE, col.justify = ALIGN, header.col.justify = ALIGN)
    addNewLine(rtffile)
    
  }
  
  done(rtffile)
  
}

SOAP_DATA <- function (CNO, DATE1, DATE2, var = c(2, 3, 5:7, 10:12), restructure = TRUE) {
  
  require(RCurl)
  
  h = basicTextGatherer()
  
  body1 = '<?xml version="1.0" encoding="utf-8"?>
  <soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
  <soap:Body>
  <Query_FEXREPORT_Data_by_LASTUPDTESTDATE_From_LisDB xmlns="http://tempuri.org/">
  <CNO>'
  
  body2 = '</CNO>
  <yyyyMMdd1>'
  body3 = '</yyyyMMdd1>\
  <yyyyMMdd2>'
  body4 = '</yyyyMMdd2>\
  <user>hd</user>\
  <passwd>xup6fup</passwd>
  </Query_FEXREPORT_Data_by_LASTUPDTESTDATE_From_LisDB>
  </soap:Body>
  </soap:Envelope>\n'
  
  body = paste0(body1,CNO,body2,DATE1,body3,DATE2,body4)
  
  tryCatch({
    curlPerform(url = "http://10.200.1.79/TsghSQL2008_ASE_Outsourcing/TSOleDB.asmx?op=Query_FEXREPORT_Data_by_LASTUPDTESTDATE_From_LisDB",
                httpheader=c(Accept="text/xml", Accept="multipart/*",
                             SOAPAction='"http://tempuri.org/Query_FEXREPORT_Data_by_LASTUPDTESTDATE_From_LisDB"',
                             'Content-Type' = "text/xml; charset=utf-8"),
                .opts = list(maxfilesize = 1e8, maxfilesize.large = 1e8, timeout = 30),
                postfields=body,
                writefunction = h$update,
                verbose = FALSE
    )},
    warning = function(msg) {},
    error = function(msg) {})
  
  if (!restructure) {return(h$value())} else {
    txt2 = h$value()
    var.names = c("TESTORDER", "NAME", "RESULT", "Unit",
                  "COLLECTIONDATE","LASTUPDTESTDATE", "LOCNAME",
                  "RESTYPE", "CHAPID", "ACCESSNUMBER",
                  "HOSTCODE", "TESTCODE")
    
    
    pos.list = list()
    for (i in var) {
      pos.list[[i*2-1]] = as.numeric(gregexpr(paste0("<",var.names[i],">"), txt2)[[1]]+2+nchar(var.names[i]))
      pos.list[[i*2]] = as.numeric(gregexpr(paste0("</",var.names[i],">"), txt2)[[1]]-1)
    }
    
    n = length(pos.list[[i*2]])
    
    if (pos.list[[i*2]][1] == -1) {return(NULL)} else {
      data = data.frame(
        報告順序 = rep(NA, n),
        報告項目簡稱 = rep(NA, n),
        結果值 = rep(NA, n),
        單位 = rep(NA, n),
        收件日期 = rep(NA, n),
        最後發報告日期 = rep(NA, n),
        開單地點 = rep(NA, n),
        檢驗結果類別 = rep(NA, n),
        檢驗類別碼 = rep(NA, n),
        開單流水號 = rep(NA, n),
        院內醫令碼 = rep(NA, n),
        檢驗碼 = rep(NA, n)
      )
      
      if (n>1) {
        
        for (i in var) {
          for (k in 1:n) {
            data[k,i] = substr(txt2, pos.list[[i*2-1]][k], pos.list[[i*2]][k])
          }
        }
        
        data = data[!is.na(data[,3]),]
        data = data[!grepl('\n', data[,3], fixed = TRUE),]
        data = data[!grepl(',', data[,3], fixed = TRUE),]
        data = data[nchar(data[,3]) <= 20,]
        data[,'收件日期'] = as.POSIXct(data[,'收件日期'])
        
        data = data[order(data[,'收件日期']),]
        
      } else {
        
        data = data[0,]
        
      }
      
      return(data[,var])
      
    }
  }
}