classdef myweibull_fit
  methods(Static)
    out=Jac(p0,p1,lambda,x,y,zz);
    out=Jac0(p0,p1,lambda,x,y,zz);
    out=Hess(p0,p1,lambda,x,y,zz);
    out=Hess0(p0,p1,lambda,x,y,zz);
        
    function R=weibull_fit(x,y,n,theta_lo,theta_hi)
      
      % [bhat, logl,xi,yi]=weibull_fit(x,y)
      % Fits the weibull model to (x,y) data set and returns the estimates bhat
      % and the log-likeliihood (logl) and interpolated xi and yi from fit
      %
      % The weibull function is P(correct)=1-0.5 exp[-(x/alpha)^beta]
      % here bhat(1)=alpha amd bhat(2)=beta
      
      options = optimoptions('fmincon','GradObj','on');
      
      xth=0;
      theta = [0.5 0.5 0.5];
      
      [X,nlogl] = fmincon(@(theta) myweibull_fit.weibull_model(theta,x,y,n,xth),theta,[],[],[],[],theta_lo,theta_hi,[],options);
      
      R.bhat=X;
      R.logl=-nlogl;
      
      R.xi=linspace(-0.8,1,100);
      % R.xi=linspace(min(x),max(x),100);
      
      R.yi=myweibull_fit.weibull_pred(R.bhat,R.xi,xth);
      
      R.xfit=x;
      R.yfit=myweibull_fit.weibull_pred(R.bhat,x,xth);
      
      R.hess=  myweibull_fit.weibull_hessian(R.bhat,x,y,n,xth);
      R.theta_lo=theta_lo;
      R.theta_hi=theta_hi;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [nlogl,dlogl]=  weibull_model(b,x,y,n,xth)
      %returns the log likelihood and the gradient
      
      p0=b(1);
      p1=b(2);
      lambda=b(3);
      
      p= myweibull_fit.weibull_pred(b,x,xth);
      
      nlogl=-n*sum(cross_surprise(y,p)+cross_surprise(1-y,1-p));
      
      %derivatives of the log-likelihood
      s=x>=xth;
      
      J(s,:)=n*myweibull_fit.Jac(p0,p1,lambda,x(s),y(s),0*y(s));
      J(~s,:)=n*myweibull_fit.Jac0(p0,p1,lambda,x(~s),y(~s),0*y(~s));
      dlogl=sum(J);
      nlogl;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function hess=  weibull_hessian(b,x,y,n,xth)
      p0=b(1);
      p1=b(2);
      lambda=b(3);
      
      s=x>=xth;
      H(s,:) =n*myweibull_fit.Hess(p0,p1,lambda,x(s),y(s),0*y(s));
      H(~s,:)=n*myweibull_fit.Hess0(p0,p1,lambda,x(~s),y(~s),0*y(~s));
      
      hess=sum(H);
      hess=reshape(hess,3,3);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function p=weibull_pred(b,x,xth)
      
      p0=1/(1+exp(b(1)));
      p1=1/(1+exp(b(2)));
      
      lambda=b(3);
      
      s=x>=xth;
      p=zeros(size(x));
      p(s) = p1+(p0-p1)*(exp(-((x(s)/lambda))));
      p(~s)= p0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    function analytic
      
      syms p0 p1 lambda x y zz real
      ep1=1/(1+exp(p1));% this ensures ep1 and ep2 are positive
      ep0=1/(1+exp(p0))
      
      p = ep1+(ep0-ep1)*(exp(-((x/lambda))));
      
      %for curve
      nlogl  = -( y * log(p) +   (1-y) * log(1-p));
      
      % for flat segment
      nlogl0 = -( y * log(ep0) +  (1-y) * log(1-ep0));
      vars=[p0 p1 lambda];
      
      J= simplify(jacobian(nlogl,  vars))
      J0=simplify(jacobian(nlogl0, vars))
      
      H= simplify(hessian(nlogl,  vars))
      H0=simplify(hessian(nlogl0,  vars))
      
      for w=1:numel(J)
        if isequaln(J(w),sym(0)), J(w)=zz; end
        if isequaln(J0(w),sym(0)), J0(w)=zz; end
      end
      
      for w=1:numel(H)
        if isequaln(H(w),sym(0)), H(w)=zz; end
        if isequaln(H0(w),sym(0)), H0(w)=zz; end
      end
      
      vars=[p0 p1 lambda x y zz]
      matlabFunction(J,'File','Jac','Vars',vars)
      matlabFunction(H(:)','File','Hess','Vars',vars)
      matlabFunction(J0,'File','Jac0','Vars',vars)
      matlabFunction(H0(:)','File','Hess0','Vars',vars)
     
    end   
  end
end
%%



