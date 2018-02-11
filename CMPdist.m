function Y=CMPdist(X,lamb,nu)
  % Conway-Maxwell-Poisson distribution

  % Y=@(X,lamb,nu) lamb.^X./(factorial(X).^z);

  summax = 100;
  termlim = 1e-6;
  sum1 = 0;
  for js = 1:summax
      term = exp(log(lamb^(js-1)) - nu*gammaln(js));
      if(js > 3)
          if((term/sum1) < termlim)
              break
          end
      end
      sum1 = sum1 + term;
  end

  %Y = lamb^X / ((factorial(n))^nu * sum);
  Y = exp(X*log(lamb) - (nu*gammaln(X+1) + log(sum1)));
