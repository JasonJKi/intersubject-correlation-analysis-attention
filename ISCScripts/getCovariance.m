% Get All Covariance.
clear all;
homeDir = '/home/jason/Repository/AttentionModulation/';
run([homeDir 'Codes/eegInfo.m'])
nModule = [2];

processIndx = 2 % choose manual or automatic
processType = eeginfo.preprocess{processIndx}
saveDir = 'Data/ISCValues/'
saveDir = [homeDir saveDir];
outFile =  [processType(1:end-1) '_cov_constrained_v2.mat'];


i=1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule};
    
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [homeDir 'Data/' eeginfo.eegProcDir processType moduletype '/' stim '.mat'];
        load(eeginDir)
        eeg(:,:,module.subjRmv{iModule})=[];
        if i == 1; eeg(:,:,end) =[];end
        [Rxy Rpool Rxy_all Rpool_all] = generate_cov(eeg,1);
        R.xy_all{i} = Rxy_all;
        R.pool_all{i} = Rpool_all;
        R.xy{i} = Rxy;
        R.pool{i} = Rpool;
        i = i+1
    end

end
save([saveDir outFile], 'R')



CovFile =  ['RPCA_allCov.mat'];
saveDir = 'ISCScripts/ISCValues/'
saveDir = [homeDir saveDir];
load([saveDir CovFile]);
nComp = 5;
% headprojection
group = {[1:4] [5 6 7 8] [1 2 5 6] [3 4 7 8] [9 10 11 12 13 14]};
labels = {'Attend', 'Disattend','Audio-Visual', 'Audio', 'Attend Audio-Visual', 'Disattend Audio-Visual', 'Attend Audio', 'Disattend Audio'};
clear A ISC

i = 1
for iGroup = 1
    RxyGrouped = zeros(64);
    RpoolGrouped = zeros(64);
    for iStim = group{iGroup}
    RxyGrouped = RxyGrouped + R.xy{iStim};
    RpoolGrouped = RpoolGrouped + R.pool{iStim};
    end
    [w(:,:,i) A(:,:,i)] = correlated_components(RxyGrouped, RpoolGrouped, .5, 64, 10);
    i = i+1;
end

headprojection(A,10,'hor',0)


