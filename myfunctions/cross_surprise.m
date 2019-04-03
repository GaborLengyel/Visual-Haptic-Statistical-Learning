function out=cross_surprise(resp,prob)
  out=zeros(size(prob));

  s=(resp>0);
  out(s)=resp(s).*log(prob(s));
