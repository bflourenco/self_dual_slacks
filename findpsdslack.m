function [Y,found]=findpsdslack(suppmatrix,dim)

% Input:
%
% suppmatrix - A 0/1 matrix that encodes a symmetric support that is the
% candidate for being the support of the slack matrix of a selfdual cone
%
% dim - the expected dimension of said cone
%
% Output:
%
% Y - an approximate psd matrix with the same dupport as suppmatrix and
% numerical rank dim
%
% found - a Boolean variable that is true if and only if Y has the required
% properties.

%Parameters:
max_tries=100;   %Maximum number of distinct random weights tried
max_iter=10000;  %Maximum number of alternate projection steps to refine the solution found
tol_eig=10^-12;  %Maximum value accepted for the dim+1 largest eigenvalue
tol_zero=10^-4;  %Minimum value accepted for the nonzero coordinates

%Initialization
found=false;
n=size(suppmatrix,1);

%Setting up the sdp constraints
X=sdpvar(n);
G=[X>=0];
for i=1:n    
    G=[G,X(i,i)==1];   
    for j=i+1:n  
        if suppmatrix(i,j)==0   
            G=[G,X(i,j)==0];
        end
    end
end

% Try for up to max_tries attempts to create a psd matrix with the given
% properties. For each attempt we generate an objective value by simply
% randomly weighting the entries of the matrix with uniformly distributed
% scalars chosen in the interval [1,2]

tries=0;
while tries<max_tries&&~found
    %Solve the optimization problem
    optimize(G,-sum(sum((1+rand(n,n)).*X)))
    Y=double(X);
    
    % If it looks like a good candidate then perform alternate projection
    % to improve the accuracy
    e=sort(eig(Y),'descend');
    iter=0;
    
    if min(min(Y+ones(n)-suppmatrix))>tol_zero
        while abs(e(dim+1))>tol_eig && iter<max_iter
            [U,S]=svd(Y);
            Y=U(:,1:dim)*S(1:dim,1:dim)*(U(:,1:dim))';
            Y=Y.*suppmatrix;
            e=sort(eig(Y),'descend');
            iter=iter+1;
        end
        
        if min(min(Y+ones(n)-suppmatrix))>tol_zero
            found=true;
        end
        tries=tries+1;
    end
end
