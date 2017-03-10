clear all;
homedir = '/home/jason/AttentionModulation/';
run([homedir 'eegInfo.m'])
nModule = [1 2 3];
processtype = eeginfo.preprocess{1}

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
            [eegClean nSampleRmv_] = preprocessing(eeg_, name, stim, eogIndx,eeginfo.fs);
            eeg(:,:,iSubj) = eegClean;
            eog(:,:,iSubj) = eeg_(:,eogIndx);
            nrmv =  nSampleRmv_/prod(size(eeg));
            nSampleRmv(iSubj) = nrmv;
            disp(['% sample removed = ' num2str(nrmv)])
        end        
        
        save([eegoutDir stim '.mat'], 'eeg','eog','nSampleRmv')
        clear eeg eog nSampleRmv
    end
end
