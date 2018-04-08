nvsd <-
function(X, y, fold=10, step.size=0.01, stop.alpha=0.05, stop.var.count=20, max.model.var.count=10, roughening.method="DCOL", do.plot=F, pred.method="MARS")
{
    
    nonlin.pred<-function(x, y, x.new, method=c("MARS", "RF", "SVM"))
    {
        
        if(is.null(nrow(x)))
        {
            x<-matrix(x, nrow=1)
            x.new<-matrix(x.new, nrow=1)
        }
        
        this.x<-as.data.frame(t(x))
        this.test.x<-as.data.frame(t(x.new))
        
        if(method=="MARS") e<-earth(x=this.x, y=y)
        else if(method=="RF")  e<-randomForest(this.x, y, ntree=max(500, 4*row(this.x)))
        else if(method=="SVM") e<-svm(x=this.x, y=y)
        ee<-predict(e, this.test.x)
        ee
    }

    find.cut<-function(x,y)
    {
        pred <- prediction(x, y)
        perf <- performance(pred,"tpr","fpr")
        opt.cut = function(perf, pred){
            cut.ind = mapply(FUN=function(x, y, p){
                d = (x - 0)^2 + (y-1)^2
                ind = which(d == min(d))[1]
                c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
                cutoff = p[[ind]])
            }, perf@x.values, perf@y.values, pred@cutoffs)
        }
        opt.cut(perf,pred)[3,1]
    }
    
    grp<-rep(1:fold, length(y))[1:length(y)]
    grp<-sample(grp, length(y), replace=FALSE)
    y.pred<-new("list")
    for(i in 1:stop.var.count) y.pred[[i]]<-matrix(NA, nrow=1, ncol=length(y))
    
    b<-stage.forward(X, y, stop.alpha=stop.alpha, stop.var.count=stop.var.count, step.size=step.size, roughening.method=roughening.method, do.plot=F)
    #    pred.ssx<-b$found.pred
    #    for(i in 1:length(pred.ssx)) pred.ssx[i]<-max(b$ssx.rec[b$sel.rec==b$found.pred[i]])
    
    if(length(b$found.pred) > 0)
    {
        
        for(m in 1:fold)
        {
            this.X<-X[,grp != m]
            this.y<-y[grp != m]
            
            for(j in 1:min(length(b$found.pred), max.model.var.count))
            {
                x.train<-this.X[b$found.pred[1:j],]
                y.train<-this.y
                x.test<-X[b$found.pred[1:j], grp == m]
                
                predicted<-nonlin.pred(x.train, y.train, x.test, method=pred.method)
                
                if(length(table(y))==2)
                {
                    predicted.0<-nonlin.pred(x.train, y.train, x.train, method=pred.method)
                    predicted<-1*(predicted>find.cut(predicted.0, y.train))
                }
                
                y.pred[[j]][grp == m]<-predicted
            }
            
        }
        
        for(j in 1:min(length(b$found.pred), max.model.var.count))
        {
            y.pred[[j]][is.na(y.pred[[j]])]<-Inf
        }
        
        ssr<-matrix(NA, nrow=1, ncol=stop.var.count)
        for(j in 1:stop.var.count)
        {
            ssr[1,j]<-sum((y-y.pred[[j]][1,])^2)
        }
        
        
        n.var.err<-apply(ssr,2,min)
        n.var.choice<-which(n.var.err == min(n.var.err,na.rm=T))[1]
        
        return(list(selected.pred=b$found.pred[1:n.var.choice], all.pred=b$found.pred))
    }else{
        
        return(list(selected.pred=NULL, all.pred=NULL))
    }
}
