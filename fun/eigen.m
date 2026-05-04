function [U,l] = eigen(X);
% usage: [U,l]=eigen(X)
% Compute the Eigenvalues and Eigenvectors of a
% semi positive definite matrix X.
% U is the matrix of the eigenvectors.
% l is the vector of the eigenvalues.
% Eigenvectors & eigenvalues are sorted in decreasing order.
% The eigenvectors are normalized: U'* U = I.
% Eigenvalues smaller than epsilon=.000001 and
% negative eigenvalues (due to rounding errors) are set to zero.
% Herve' Abdi, September 1990.
    epsilon=.000001;
%  tolerance to be considered 0 for an eigenvalue
   [U,D]=eig(X);
   
%    if isreal(U) ==0
%        [U,D]=eigs(X);
%    end
   D=diag(D);
   [l,k]=sort(D);
   n=length(k);
    l=l((n+1)-(1:n));
    U=U(:,k((n+1)-(1:n)));
% keep the non-zero eigen value only (tolerance=epsilon)
  pos=find(any([l';l'] > epsilon ));
  l=l(pos);
  U=U(1:n,pos);
% Normalize U
  U=U./( ones(n,1) * sqrt(sum(U.^2) ) )  ;

