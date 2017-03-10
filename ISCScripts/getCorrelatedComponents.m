clear all;
homeDir = '/home/jason/Repository/AttentionModulation/';
run([homeDir 'Codes/eegInfo.m'])
nModule = [1 2 3];
saveDir= [homeDir 'Data/ISCValues/'];
load([saveDir  'RPCA_cov_constrained_v2.mat'], 'R')

groupIndx{1} = [1 2 3 4 5 6 7 8]
groupIndx{2} = [9 10 11 12 13 14]
nGroup = 2;

%correlated component constants
D = 64; gamma =.5; whitening = 0

for iGroup = 1:nGroup;
    Rxy = zeros(64); Rpool=zeros(64);
    for iStim = groupIndx{iGroup}
        Rxy = Rxy + R.xy{iStim};
        Rpool = Rpool + R.pool{iStim};
    end
    [W{iGroup} A{iGroup}] = correlated_components(Rxy, Rpool, gamma, D, whitening);
end

save([homeDir saveDir 'allComponents.mat'], 'W', 'A')

