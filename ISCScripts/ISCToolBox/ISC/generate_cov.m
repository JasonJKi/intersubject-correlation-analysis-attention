function [Rxy Rpool RxyAll RpoolAll] = generate_cov(x,average)
    % [Rxy Rpool] = generate_cov(x,average) %%
    % x = samplex x Channels x subjects
    % if average > 0 compute mean covariance
    subjkept = 1:size(x,3);
    [T,D,Nsubj]=size(x);
    Rpool = zeros(D,D,Nsubj);
    RpoolAll = zeros(D);
    for ii = 1:Nsubj;
    Rpool(:,:,ii) = cov(x(:,:,ii));
    RpoolAll = RpoolAll + Rpool(:,:,ii);
    end
    
    subjects = 1:length(subjkept);
    permutationList = combnk(subjects,2);
    
    % this line makes Rxx==Ryy so that the optimization algorithm is exact
    Nperm = length(permutationList);
    
    Rxy = zeros(D*2);
    RxyAll = zeros(D);
    for ii = 1:Nperm;
        Rxy(:,:,ii) = cov([x(:,:,permutationList(ii,1)), x(:,:,permutationList(ii,2))]);
        RxyAll = RxyAll + Rxy(1:D,D+1:2*D,ii);
    end
    
    if average
        RxyAvg = RxyAll;
        RxyAvg = (RxyAvg + RxyAvg');
        RxyAvg = RxyAvg ./ (2*Nsubj*(Nsubj-1));
        Rxy = RxyAvg;
        RpoolAvg = RpoolAll ./ (Nsubj*2);
        Rpool = RpoolAvg;
    end
end