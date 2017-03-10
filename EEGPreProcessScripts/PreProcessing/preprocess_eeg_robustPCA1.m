clear all
run /home/jason/Experiments/Backwardcount/codes/bdf_reader/eegdetail.m
nameCov = ['cov_almLamda/']
fsref = 256; % sampling rate
sRemove = 3; % number standardation removal
nRemove = 4; % iteration of removal
dim = 64;    % # channels
nLamda = 5   % # of lamda to filter eeg
unProcessedLoadPath ='../../EEG/UnProcessed/';
processedSavePath = '../../EEG/Processed/RPCAProcessed_v5/';
comment = []
experiment = {    'Freeview'    'Fixation'    'Audio'    };
gamma = .5;grassman = 0; Ncomp = 15

iFigure = 1
for iFileset = [1 3];
    fileSetLoaPath = [unProcessedLoadPath experiment{iFileset} '/'];
    fileSetSavePath = [processedSavePath experiment{iFileset} '/'];
    mkdir(fileSetSavePath)
    
    stims=Stims{iFileset};stimsave=Stimsave{iFileset};
    subjname=Subjname{iFileset};eogchannel=Eogchannel{iFileset};
    
    nStim = length(stims);
    EEG = []
    for iStim = 1:nStim;        
        eegclean=[];eeg=[];eog=[];eegunprocessed=[];eegreg =[];
        outlier=[];Rxy=[];Rpool=[];

        nSubj = Nsubj{iFileset}(iStim);
        i = 1
        for iSubj =1:nSubj;
            name = subjname{iSubj};
            disp([name '_' num2str(iSubj) ' ' stims{iStim}])
            x = load([fileSetLoaPath name '.mat'], stims{iStim});
            eval(['eeg' ['=x.' stims{iStim} ';'] ]);
            eegunprocessed = eeg(:,1:64);
            eog(:,:,i) = eeg(:,eogchannel); % unprocessed eog
            eegreg = eegunprocessed-eog(:,:,i)*(eog(:,:,i)\eegunprocessed); %eog regression
            EEG(:,:,i) = inexact_alm_rpca(eegreg);
            i = i+1;
        end
        
            save([fileSetSavePath stimsave{iStim} comment], 'eeg')
            clear eeg
    end
    clear nSubj stimlength
end 

load('jpn_free.mat')
x= eeg;
x(:,:,end-1:end) = [];
load('jpn_freeadd3.mat')
eeg= cat(3,x,eeg);
save(['../../RPCAProcessed_v4/Freeview/' 'jpn_free'],'eeg')

