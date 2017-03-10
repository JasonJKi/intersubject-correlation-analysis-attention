function [Rxy Rpool] = covariance_leave_one_out(Rxy_all,Rpool_all,leaveout,Nsubj);
%[Rxy Rpool] = covariance_leave_one_out(Rxy_all,Rpool_all,leaveout);
[D,~,~]=size(Rpool_all);
subjects = 1:Nsubj; %
permList = combnk(subjects,2);
Nperm = size(Rxy_all,3);
[i j] = find((permList - (permList == leaveout)*leaveout) == 0);
covlist = 1:Nperm;
covlist(i) = [];
N = length(covlist);

%leaveoneout method.
Rxy = zeros(D);
for ii = 1:N;
    Rxy = Rxy + Rxy_all(:,:,covlist(ii));
    Rpool = Rpool_all(:,:,covlist(ii));
end
Rxy = Rxy + Rxy';
