% compute the ISC.
clear all;
homeDir = '/home/jason/Repository/AttentionModulation/';
run([homeDir 'Codes/eegInfo.m'])
nModule = [1 2 3];
processtype = eeginfo.preprocess{2}
currentDir = 'Data/ISCValues/'
load([homeDir currentDir processtype(1:end-1) '_allComponents.mat'],'W')
compIndx = [1 1 2]; nComp = 10;

i = 1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule}
    
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [homeDir 'Data/' eeginfo.eegProcDir processtype moduletype '/' stim '.mat'];
        load(eeginDir)
        eeg(:,:,module.subjRmv{iModule}) = [];
        nSubj= size(eeg,3);
        for iSubj = 1:nSubj
        isc(:,iSubj) = concat_matrix_ISC(eeg,eeg,W{compIndx(iModule)},iSubj,nComp);
        disp(iSubj)
        end
        ISC{i} = isc;
        i= i+1;
        clear isc
    end
end
save([homeDir currentDir processtype(1:end-1) '_allISC_v2.mat'], 'ISC')