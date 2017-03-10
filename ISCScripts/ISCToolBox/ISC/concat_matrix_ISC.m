function [isc ccx ccy] = concat_matrix_ISC(x,y,w,subj,Ncomp,filter)
% concat eeg matrix, project, and compute isc 
% isc = concat_matrix_ISC(x,y,w,subj,Ncomp)

if nargin < 6; filterStatus=0;else filterStatus=1;end

[T,D,nsubj]=size(x);
permutationList = [ones(nsubj,1)*subj, [1:nsubj]'];
permutationList(subj,:) = [];
Nperm = size(permutationList,1);

T = size(x, 1);
X = zeros(size(x,1),D,Nperm);
Y = zeros(size(y,1),D,Nperm);

for i = 1:Nperm
    X(:,:,i) = x(:,:,permutationList(i,1));
    Y(:,:,i) = y(:,:,permutationList(i,2));
end
X = reshape(permute(X,[1 3 2]),[T*Nperm D]);
Y = reshape(permute(Y,[1 3 2]),[T*Nperm D]);

% Finally the Correlation

% correlated component time courses
ccx = X*w(:,1:Ncomp);
ccy = Y*w(:,1:Ncomp);

if filterStatus
    ccx = filtfilt(filter.b,filter.a,ccx);
    ccy = filtfilt(filter.b,filter.a,ccy);
end

% finally the ISC!
r = corrcoef([ccx ccy]);
r = r(1:Ncomp,Ncomp+1:end);
isc = diag(r);
end