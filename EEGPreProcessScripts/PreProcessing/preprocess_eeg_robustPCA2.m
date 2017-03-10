clear all;
homedir = '/home/jason/AttentionModulation/';
run([homedir 'eegInfo.m'])
nModule = [1 2];
processtype = 'RPCA2/'

for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule}
    
    for iStim = 1:nStim;
        nSubj = module.nSubj{iModule}(iStim);
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [eeginfo.eegUnprocDir moduletype '/'];
        eegoutDir = [homedir eeginfo.eegProcDir processtype moduletype '/' ];
        mkdir(eegoutDir);
        
        for iSubj =1:nSubj;
            name = ['subject' num2str(iSubj) '.mat'];
            disp([stim name])
            load([eeginDir name], stim);
            eval(['eeg_' ['=' stim ';'] ]);
            eog(:,:,iSubj) = eeg_(:,eogIndx);
            eegreg = eeg_(:,1:64)-eog(:,:,iSubj)*(eog(:,:,iSubj)\eeg_(:,1:64)); %eog regression
            [eegClean nSampleRmv_] = preprocessing(eeg_, name, stim, eogIndx,eeginfo.fs);
            eegClean_1 = inexact_alm_rpca(eegreg);
            eegClean_2 = inexact_alm_rpca(eegClean);
            [eegClean nSampleRmv_] = preprocessing(eeg_, name, stim, eogIndx,eeginfo.fs);
            eeg(:,:,iSubj) = eegClean;
            figure(1) 
        end        
        save([eegoutDir stim '.mat'], 'eeg','eog')
        clear eeg eog
    end
end
