%{ sample code from the 2012 Dmochowski Paper. Revised and condensed for learning.
% a. preprocesing
% 1. drift and 60hz filter
% 2. eye artifact removal regression
% 3. artifact rejection
% 
% check
% - spectral validation
% - threshold maps
% 
% b. correlated compoent analysis
% 1. compute pool and cross covariance
% 2. regularization of the covariance
% 3. InterSubject Correlation
% 
% check
% - Forward model
%}

Nremove = 4; % number of iterations for artefact removal
Sremove = 3; % number of std for artefact removal

%% Preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% filtering
fs=256; %sampling rate
% generate filter
[hpnum,hpdenom]=butter(5,0.5/fs*2,'high'); % drift removal
[notchnum,notchdenom]=butter(5,[59 61]/fs*2,'stop'); % 60Hz line noise
a = poly([roots(hpdenom);roots(notchdenom)]);
b = conv(hpnum,notchnum);
eeg=filter(b,a,eeg); % HP filter for drift removal

%% regress out eog activity (eye movement activity)
eeg = eeg - eog*(eog\eeg);
% show the "clean" eeg
subplot(2,2,2); imagesc(eeg'); title; caxis([-50 50])

%% outliers rejection
%Coming for bandpassed Data, TBD)
mask = zeros(size(eeg));
eeg = eeg';
tempeeg = removeOutliers(eeg, 3, 3, round(fs*.050)); %
mask =(eeg==tempeeg);
maskHolder{i,1} = mask;
eeg = tempeeg';
imagesc(mask>0); title('artifacts')

% specral validation
for chan=1:size(EEG,2)
    subplot(3,2,5)
    psd(EEG(:,chan,i),fs,fs)
    title([files.name ', chan=' num2str(chan)])
    drawnow
    pause
end

%% Correlated Compoents Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute cross and pooled covariance
% construct all possible pairs
permutationList = combnk(1:Nsubj,2);
Nperm = size(permutationList,1);

% memory-efficient computation of the  pooled covariance matrix,
% avoiding concatenation of all pairs of datasets
Rpool = zeros(D);
for ii = 1:Nsubj
    Rpool = Rpool + cov(EEG(:,:,ii));
end
Rpool = Rpool ./ Nsubj;

% memory-efficient computation of the cross-covariance matrix,
% avoiding concatenation of all pairs of datasets
Rxy = zeros(D);
for i = Nperm:-1:1
    Rxy_ = cov([EEG(:,:,permutationList(i,1)), EEG(:,:,permutationList(i,2))]);
    Rxy = Rxy + Rxy_(1:D,D+1:2*D);
end
Rxy = Rxy + Rxy';
Rxy = Rxy ./ (Nsubj*(Nsubj-1));

%% regularization
if gamma
    % shrinkage
    Rpool = (1-gamma)*Rpool + gamma*mean(eig(Rpool))*eye(D);
    [W,L] = eig(Rxy, Rpool);   [d,indx]=sort(diag(L),'descend'); W = W(:,indx);
else
    % regularization of pooled covariance
    K = 10; % how many PCs to keep
    [V,L] = eig(Rpool); [d,indx]=sort(diag(L),'descend'); V = V(:,indx);
    Rxy=V(:,1:K)*diag(1./d(1:K))*V(:,1:K)'*Rxy;
    [W,L] = eig(Rxy);   [d,indx]=sort(diag(L),'descend'); W = W(:,indx);
end
W = W(:,1:Ncomp);

%% InterSubject Correlation
permutationList = [permutationList; permutationList(:,[2 1])]; 
Nperm = size(permutationList,1);
for i = Nperm:-1:1
    x(:,:,i) = EEG(:,:,permutationList(i,1));
    y(:,:,i) = EEG(:,:,permutationList(i,2));
end
x = reshape(permute(x,[1 3 2]),[(T)*Nperm D]);
y = reshape(permute(y,[1 3 2]),[(T)*Nperm D]);

% correlated component time courses
ccx = x*W; 
ccy = y*W; 
% finally the ISC!
[r,p] = corrcoef([ccx ccy]); r=r(1:Ncomp,Ncomp+1:end);
ISC = diag(r)

% check
% forward model
A=Rpool*W*inv(W'*Rpool*W);
maplimits=[min(A(:)) max(A(:))];
figure;
for k=1:Ncomp
    subplot(1,Ncomp,k); topoplot( A(:,k)  ,'BioSemi32.loc','numcontour',5,'plotrad',0.5,'electrodes','off','maplimits',maplimits);
end


