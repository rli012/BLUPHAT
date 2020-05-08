
blup.bayes<-function(x,y,kk){
   loglike<-function(parm){
      v<-diag(n)*parm[p+1]
      vv<-0
      for(k in 1:p){
         vv<-vv+parm[k]
         v<-v+kk[[k]]*vv
      }
      vi<-solve(v)
      xx<-t(x)%*%vi%*%x
      xy<-t(x)%*%vi%*%y
      b<-solve(xx,xy)
      d1<-unlist(determinant(v))
      d1<-d1[1]
      d2<-unlist(determinant(xx))
      d2<-d2[1]
      r<-y-x%*%b
      q<-t(r)%*%vi%*%r
      loglike<- -0.5*(d1+d2+q)
      return(-loglike)
      
   }
   
   fixed<-function(parm){
      v<-diag(n)*parm[p+1]
      g<-matrix(0,n,n)
      vv<-0
      for(k in 1:p){
         vv<-vv+parm[k]
         g<-g+kk[[k]]*vv
      }
      v<-v+g
      vi<-solve(v)
      xx<-t(x)%*%vi%*%x
      xy<-t(x)%*%vi%*%y
      covb<-solve(xx)
      beta<-solve(xx,xy)
      yhat<-g%*%vi%*%(y-x%*%beta)
      yobs<-y-x%*%beta
      r2<-cor(yobs,yhat)^2
      result<-list(beta,covb,r2)
      return(result)
   }
   
   
   loglike0<-function(x,y){
      xx<-t(x)%*%x
      xy<-t(x)%*%y
      b<-solve(xx,xy)
      r<-y-x%*%b
      s2<-drop(t(r)%*%(r))/(n-ncol(x))
      v<-diag(n)*s2
      vi<-diag(n)/s2
      d1<-unlist(determinant(v))
      d1<-d1[1]
      d2<-unlist(determinant(xx))
      d2<-d2[1]
      q<-t(r)%*%vi%*%r
      loglike<- -0.5*(d1+d2+q)
      return(-loglike)
   }


   
   n<-length(y)
   p<-length(kk)
   
   fn0<-loglike0(x,y)
   
   parm0<-rep(1,p+1)
   result<-optim(par=parm0,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=1e-5,upper=1e5)
   parm<-result$par
   conv<-result$convergence
   fn<-result$value
   lrt<-2*(fn0-fn)
   hess<-result$hessian
   covp<- solve(hess)
   bb<-fixed(parm)
   beta<-bb[[1]]
   covb<-bb[[2]]
   r2<-bb[[3]]
   fixed<-data.frame(conv,fn0,fn,lrt,beta,covb,r2)
   v1<-parm[1]
   v2<-parm[2]
   ve<-parm[3]
   parm<-data.frame(v1,v2,ve)
   
   result<-list(fixed,parm)
   return(result)
   
}
