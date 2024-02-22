function [x, error, iter, flag,condn,tpre] = pcg0_cond(A, x, b, max_it, tol,pre)

tpres=tic;
  flag = 0;                                 % initialization
  iter = 0;

  bnrm2 = norm( b-A*x );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  fprintf('PCG residual(%d) = %2.2e\n',iter,error)
  if ( error < tol ) return, end

  for iter = 1:max_it                       % begin iteration
     z=pre(r);
     rho = (r'*z);

     if ( iter > 1 )                       % direction vector
        beta(iter) = rho / rho_1;
        p = z + beta(iter)*p;
     else
        p = z;
     end

     q = A*p;
     alpha(iter) = rho / (p'*q );
     x = x + alpha(iter) * p;                    % update approximation vector

     r = r - alpha(iter)*q;                      % compute residual
     error = norm( r ) / bnrm2;            % check convergence
     fprintf('PCG residual(%d) = %2.2e\n',iter,error)
     if ( error <= tol ), break, end

     rho_1 = rho;

  end
%   condn=1;
% return
%  compute eigenvalues and codition number
tpre=toc(tpres);
     d(1)= 1/alpha(1);
     for i = 2: iter
     d(i)=beta(i)/alpha(i-1)+1/alpha(i);
     end
     for i = 1: iter-1
     s(i)=-1*sqrt(beta(i+1))/alpha(i);
     end
%
     T = sparse(iter,iter);
     T(1,1)= d(1);
     for i=2:iter
     T(i,i) = d(i);
     end
     for i = 1:iter-1
     T(i,i+1)= s(i); T(i+1,i) = T(i,i+1);
     end
     lambda = eig(full(T));
     lambdamax = max(lambda);
     lambdamin = min(lambda);
     fprintf('Condition number and iteration number are %2.2e %d\n', lambdamax/lambdamin,iter);
     condn=lambdamax/lambdamin;
  if ( error > tol ) flag = 1; end         % no convergence

% END cg.m