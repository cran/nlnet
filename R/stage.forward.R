stage.forward <-
function(X, y, step.size=0.01, stop.alpha=0.01, stop.var.count=20, roughening.method="DCOL", tol=1e-8, spline.df=5, dcol.sel.only=FALSE, do.plot=F)
{
    
    find.nonlin.assoc<-function(x,y)
    {
        o<-order(x)
        x<-x[o]
        y<-y[o]
        
        d<-diff(y)
        s.epsilon<-var(d)/2
        
        l<-length(y)
        ssy<-sum((y-mean(y))^2)
        ssr<-s.epsilon*(l-1)
        #    if(ssr>ssy) ssr<-ssy
        ssx<-ssy-ssr
        
        #    message(c("SSY: ", ssy, "  ,SSX: ", ssx, "  ,SSR: ", ssr))
        #    message(c("% variance attributed to X: ", ssx/ssy))
        
        to.return<-c(ssx/ssy, ssy, ssx, ssr)
        names(to.return)<-c("%explained", "SSY", "SSX", "SSR")
        return(to.return)
    }
    
    null.gen<-function(y, B=1000)
    {
        l<-length(y)
        s.y<-var(y)
        
        r<-rep(0,B)
        for(i in 1:B)
        {
            o<-sample(l)
            this.y<-y[o]
            this.d<-diff(this.y)
            r[i]<-var(this.d)/2
        }
        
        (s.y-r)/s.y
    }

    find.assoc.batch.o<-function(o,y, B=1000)  # o is a matrix. Every row is orders of an X variable
    {
        l<-length(y)
        ssy<-sum((y-mean(y))^2)
        y0<-y
        
        ssr<-rep(0, nrow(o))
        
        for(i in 1:nrow(o))
        {
            y<-y0[o[i,]]
            d<-diff(y)
            ssr[i]<-var(d)/2*(l-1)
        }
        
        ssx<-ssy-ssr
        explained<-ssx/ssy
        
        r<-null.gen(y, B=B)
        m<-mean(r)
        s<-sd(r)
        p.nonlin<-pnorm(explained, mean=m, sd=s, lower.tail=FALSE)
        
        return(cbind(p.nonlin, explained, rep(ssy, length(p.nonlin)), ssx, ssr))
    }

    find.assoc.batch.l<-function(X, y)
    {
        n<-length(y)
        r<-cor(t(X), y)
        ttt<-r*sqrt((n-2)/(1-r^2))
        ttt[is.na(ttt)]<-Inf
        
        p<-ttt
        p[ttt>=0]<-pt(ttt[ttt>=0], df=n-2, lower.tail=FALSE)
        p[ttt<0]<-pt(ttt[ttt<0], df=n-2, lower.tail=TRUE)
        p*2
    }

    one.step.spline<-function(x,y, r=0.1, spline.df, do.plot=F)
    {
        if(do.plot) plot(x,y)
        l<-length(y)
        
        smooth.tol<-1e-6 * IQR(x)
        if(smooth.tol==0) smooth.tol<-1e-6 * diff(range(x))
        f<-smooth.spline(x,y, df=min(l/10, spline.df), tol=smooth.tol)
        if(do.plot) lines(f)
        
        fitted<-predict(f,x)$y
        
        delta<-fitted-y
        
        #y.new<-y-delta/sd(delta)*r*sd(y)
        y.new<-y-delta*r
        
        
        if(do.plot) points(x,y.new,col="red")
        
        y.new
    }
    
    one.step.o<-function(o,y, r=0.1, do.plot=F)  # o is the order based on a single X
    {
        if(do.plot) plot(o,y)
        
        y<-y[o]
        l<-length(y)
        
        y.new<-y
        y.new[2:(l-1)]<-y[2:(l-1)]*(1+2*r) - y[1:(l-2)]*r - y[3:l]*r
        y.new[1]<-y[1]*(1+r) - y[2]*r
        y.new[l]<-y[l]*(1+r) - y[l-1]*r
        
        delta<-y-y.new
        y.new<-y-delta/sd(delta)*r*sd(y)
        
        y.new<-y.new[order(o)]
        if(do.plot) points(o,y.new,col="red")
        y.new
    }
    
    
    y<-(y-mean(y))/sd(y)
    y.0<-y
    
    orders<-X
    for(i in 1:nrow(X))
    {
        orders[i,]<-order(X[i,])
    }
    
    orig.assoc.o<-find.assoc.batch.o(orders,y)
    if(dcol.sel.only)
    {
        orig.assoc.p<-orig.assoc.o[,1]
    }else{
        orig.assoc.l<-find.assoc.batch.l(X,y)
        orig.assoc.p<-cbind(orig.assoc.o[,1], orig.assoc.l)
        orig.assoc.p<-apply(orig.assoc.p, 1, min)
        orig.assoc.p<-(2-orig.assoc.p)*orig.assoc.p
    }
    
    orig.assoc<-orig.assoc.o
    orig.assoc[,1]<-orig.assoc.p
    
    assoc<-orig.assoc
    sig<-0
    p.rec<-sel.rec<-ssx.rec<-rep(0, 1e5)
    ptr<-1
    
    while(sig <= stop.alpha & length(unique(sel.rec[1:ptr])) <= stop.var.count+1)
    {
        sel<-which(assoc[,1]==min(assoc[,1]))[1]
        sig<-assoc[sel,1]
        
        y.old<-y
        if(roughening.method=="DCOL")
        {
            y<-one.step.o(orders[sel,], y, r=step.size, do.plot=do.plot)
        }else{
            y<-one.step.spline(X[sel,], y, r=step.size, spline.df=spline.df, do.plot=do.plot)
        }
        y<-(y-mean(y))/sd(y)
        
        # message(y)
        
        cat(" | ")
        cat(signif(sig,2))
        delt<-sum((y.old-y)^2)/sum(y.old^2)
        cat(" ")
        cat(signif(delt,2))
        if(delt <= tol) break
        
        this.new.assoc<-find.nonlin.assoc(X[sel,],y)
        ssx.rec[ptr]<-(orig.assoc[sel,4]-this.new.assoc[3])/orig.assoc[sel,3]
        sel.rec[ptr]<-sel
        p.rec[ptr]<-sig
        
        assoc.o<-find.assoc.batch.o(orders, y)
        if(dcol.sel.only)
        {
            assoc.p<-assoc.o[,1]
        }else{
            assoc.l<-find.assoc.batch.l(X,y)
            assoc.p<-cbind(assoc.o[,1], assoc.l)
            assoc.p<-apply(assoc.p, 1, min)
            assoc.p<-(2-assoc.p)*assoc.p
        }
        assoc<-assoc.o
        assoc[,1]<-assoc.p
        
        to.replace<-which(!(1:nrow(X) %in% sel.rec[1:ptr]))
        orig.assoc[to.replace,]<-assoc[to.replace,]
        ptr<-ptr+1
        
    }
    ptr<-max(1, ptr-2)
    #ptr<-ptr-2
    
    if(ptr>=1)
    {
        ssx.rec<-ssx.rec[1:ptr]
        sel.rec<-sel.rec[1:ptr]
        p.rec<-p.rec[1:ptr]
        
        found.x<-unique(sel.rec)
        plot(0:ptr, c(0, ssx.rec), type="n", xlab="steps")
        for(i in 1:length(found.x))
        {
            this.x<-(1:ptr)[sel.rec==found.x[i]]
            this.x<-c(min(this.x)-1, this.x)
            this.y<-c(0, ssx.rec[sel.rec==found.x[i]])
            lines(this.x, this.y)
            points(this.x, this.y,cex=.4)
        }
    }else{
        found.x<-NULL
        ssx.rec<-NULL
        sel.rec<-NULL
        p.rec<-NULL
    }
    return(list(found.pred=found.x, ssx.rec=ssx.rec, sel.rec=sel.rec, p.rec=p.rec))
}
