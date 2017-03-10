function [W A] = correlated_components(Rxy, Rpool, gamma, D, whitening)
% [W A] = correlated_components(Rxy, Rpool, gamma, D, whitening)

if whitening
    % regularization of pooled covariance
    K = whitening; % how many PCs to keep
    [V,L] = eig(Rpool); [d,indx]=sort(diag(real(L)),'descend'); V = V(:,indx);
    RpoolReg = V(:,1:K)*diag(1./d(1:K))*V(:,1:K)'
    Rxy= RpoolReg*Rxy;
    [W,L] = eig(Rxy);   [d,indx]=sort(diag(real(L)),'descend'); W = W(:,indx);
else 
    % use shrinkage method to regularize pooled covariance
    if gamma > 0
    Rpool = (1-gamma)*Rpool + gamma*mean(eig(Rpool))*eye(D);
    end
    [W,L] = eig(Rxy, Rpool); % = Rxy*inv(Rpool)
    [d,indx]=sort(diag(real(L)),'descend'); W = W(:,indx); % arrange cc by variance 
end
A=Rpool*W*inv(W'*Rpool*W); % Compute Forward Model
