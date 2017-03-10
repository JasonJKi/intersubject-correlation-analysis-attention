clear all;
unprocessedDir = '/home/jason/Repository/AttentionModulation/Data/EEG/Unprocessed/OrderReversed/';
processedDir= '/home/jason/Repository/AttentionModulation/Data/EEG/Processed/OrderReversed_v2/';
mkdir(processedDir);
nModule = [1 2];
processtype = 'RPCA/'
stimNames={'bng_1' 'bng_2';'gbu_1' 'gbu_2';'pieOrg_1' 'pieOrg_2'}

nStim = 3;
eogIndx = {[65:66 69],[65:66 69],65:70};
for iStim = 1:nStim;
    X11=[]
    for iSubj =1:20;
        name = ['Subject' num2str(iSubj) '.mat'];
        load([unprocessedDir name], stimNames{iStim,1});
        eval(['eeg_' ['=' stimNames{iStim,1} ';'] ]);
        eog = eeg_(:,eogIndx{iStim});
        eegreg = eeg_(:,1:64)-eog*(eog\eeg_(:,1:64)); %eog regression
        X11=[X11;eegreg];
    end
    X22=[]
    for iSubj =1:20;
        name = ['Subject' num2str(iSubj) '.mat'];
        load([unprocessedDir name], stimNames{iStim,2});
        eval(['eeg_' ['=' stimNames{iStim,2} ';'] ]);
        eog = eeg_(:,eogIndx{iStim});
        eegreg = eeg_(:,1:64)-eog*(eog\eeg_(:,1:64)); %eog regression
        X22=[X22;eegreg];
    end
    X=[X11;X22];
    eegClean = inexact_alm_rpca(X);
    
    X1=eegClean(1:end/2,:);
    X2=eegClean((end/2)+1:end,:);
    eegLength=length(eegreg);
    for iSubj =1:20;
        eeg(:,:,iSubj)=X1((iSubj-1)*eegLength+1:(iSubj)*eegLength,:);
        imagesc(eeg(:,:,iSubj)')
    end
    save([processedDir stimNames{iStim,1}], 'eeg')
    
    for iSubj =1:20;
        eeg(:,:,iSubj)=X2((iSubj-1)*eegLength+1:(iSubj)*eegLength,:);
    end
    save([processedDir stimNames{iStim,2}], 'eeg')
    clear eeg
end

%         
%         
%             eval(['eeg_' ['=' stim ';'] ]);
%             eog(:,:,iSubj) = eeg_(:,eogIndx);
%             eegreg = eeg_(:,1:64)-eog(:,:,iSubj)*(eog(:,:,iSubj)\eeg_(:,1:64)); %eog regression
%             [eegClean nSampleRmv_] = preprocessing(eeg_, name, stim, eogIndx,eeginfo.fs);
%             eegClean_1 = inexact_alm_rpca(eegreg);
%             eegClean_2 = inexact_alm_rpca(eegClean);
%             [eegClean nSampleRmv_] = preprocessing(eeg_, name, stim, eogIndx,eeginfo.fs);
%             eeg(:,:,iSubj) = eegClean;
%             figure(1) 
%         end        
%         save([eegoutDir stim '.mat'], 'eeg','eog')
%         clear eeg eog
%     end
% end