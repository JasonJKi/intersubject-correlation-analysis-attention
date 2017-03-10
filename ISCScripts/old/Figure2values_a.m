% Generate figure 2.
clear all

file = 1:3;
run ../../Codes/BDFRead/EEGDetail.m

experiment = { 'Freeview'    'Fixation'    'Audio' };
EEGLoadPath = '../../EEG/Processed/';
Version = 'RPCAProcessed_v4/'
EEGLoadPath = [EEGLoadPath Version]
eegname = 'eeg'
nComp = 5;nGamma = 1;grassman = 0;dim = 64; gamma = .9;
saveDir ='Figure_Final_1'; %mkdir(saveDir)
iFile = 1;
iStim = 1;
wIndx = [1 2 3];

load('cc_v1.mat')

for  Files = [1 2 3];
    stims=Stims{Files};stimsave=Stimsave{Files};
    subjname=Subjname{Files};eogchannel=Eogchannel{Files};
    attend = Attend{Files}; disAttend = Disattend{Files};
    subjrmv = Indsubjrmv{Files};
    
    nstims = length(stims);
    for stim = 1:nstims;
        % load attend viewing group
        load([EEGLoadPath experiment{Files} '/' stimsave{stim} '.mat'],eegname);eval(['x=' eegname ';' ]);
        x(:,:,subjrmv) =[]; y =[];
        % [Rxy Rpool] = generate_all_cov(x);
        Rxy = []; Rpool=[];
        nSubj = size(x,3);
        w_ = w(:,:,wIndx(Files));
        % load covariance of attend (normative) group
        isc_individual = compute_individual_subject_ISC_v2(x,y,Rpool,Rxy,nComp,gamma,dim,nSubj,w_);
        iscs{iStim} = isc_individual;
        iStim = iStim + 1
        clear isc_free isc_count
    end
    iFile = iFile+1
end


save([saveDir],'iscs')