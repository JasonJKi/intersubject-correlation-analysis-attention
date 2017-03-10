% Generate figure 2 (test).
clear all
file = 1:3;
run ../../Codes/BDFRead/EEGDetail.m
experiment = { 'Freeview'    'Fixation'    'Audio' };
EEGLoadPath = '../../EEG/Processed/';
eegname = 'eeg'
Version = 'RPCAProcessed_v3/'
EEGLoadPath = [EEGLoadPath Version]
nComp = 10;nGamma = 1;grassman = 0;dim = 64; gamma = .5;
 saveDir ='Figure_2_iscs_test_v2'; 

for iFile = file
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogchannel=Eogchannel{iFile};
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    subjrmv = Indsubjrmv{iFile};
        
    nStim = length(stims)/2;
    for iStim = 1:nStim;
        % load attend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{disAttend(iStim)} '.mat'],eegname);eval(['x=' eegname ';' ]);
        [Rxy Rpool] = generate_all_cov(x);
        % load disattend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{attend(iStim)} '.mat'],eegname);eval(['y=' eegname ';']);
        
        nSubj=size(x,3);
        % load covariance of attend (normative) group
        [isc_free isc_count] = compute_individual_subject_ISC(x,y,Rpool,Rxy,nComp,gamma,dim,nSubj);
        iscs{iFile,iStim}.free = isc_free;
        iscs{iFile,iStim}.count = isc_count;
        clear isc_free isc_count
    end   
end
save([saveDir],'iscs')