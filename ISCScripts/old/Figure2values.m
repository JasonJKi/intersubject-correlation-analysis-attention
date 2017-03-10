% Generate figure 2.
clear all

file = 1:3;
run ../../Codes/BDFRead/EEGDetail.m

experiment = { 'Freeview'    'Fixation'    'Audio' };
EEGLoadPath = '../../EEG/Processed/';
Version = 'RPCAProcessed_v4/'
EEGLoadPath = [EEGLoadPath Version]
eegname = 'eeg'
nComp = 10;nGamma = 1;grassman = 0;dim = 64; gamma = .9;
 saveDir ='Figure_2_v3'; %mkdir(saveDir)
    iFile = 1
for  Files = [1 3]

    stims=Stims{Files};stimsave=Stimsave{Files};
    subjname=Subjname{Files};eogchannel=Eogchannel{Files};
    attend = Attend{Files}; disAttend = Disattend{Files};
    subjrmv = Indsubjrmv{Files};
        
    nStim = length(stims)/2;
    for iStim = 1:nStim;
        % load attend viewing group
        load([EEGLoadPath experiment{Files} '/' stimsave{attend(iStim)} '.mat'],eegname);eval(['x=' eegname ';' ]);
        [Rxy Rpool] = generate_all_cov(x);
        % load disattend viewing group
        load([EEGLoadPath experiment{Files} '/' stimsave{disAttend(iStim)} '.mat'],eegname);eval(['y=' eegname ';']);
        
        nSubj=size(y,3);
        % load covariance of attend (normative) group
        [isc_free isc_count] = compute_individual_subject_ISC(x,y,Rpool,Rxy,nComp,gamma,dim,nSubj);
        iscs{iFile,iStim}.free = isc_free;
        iscs{iFile,iStim}.count = isc_count;
        clear isc_free isc_count
    end   
    iFile = iFile + 1
end
save([saveDir],'iscs')