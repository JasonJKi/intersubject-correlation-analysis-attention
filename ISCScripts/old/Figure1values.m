% Generate figure 1.
clear all
file = 1;
run ../../Codes/BDFRead/EEGDetail.m
experiment = { 'Freeview'    'Fixated'    'Audio' };
EEGLoadPath = '../../EEG/Processed/';
eegname = 'eeg'
Version = 'RPCAProcessed_v3/'
EEGLoadPath = [EEGLoadPath Version]
nComp = 15;nGamma = 1;grassman = 0;dim = 64; gamma = .9;
saveDir ='Figure_1_iscs_v2_freeview'; 

for iFile = file
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogchannel=Eogchannel{iFile};
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    subjrmv = Indsubjrmv{iFile};
        
    nStim = length(stims)/2;
    for iStim = 1:nStim;
        % load attend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{attend(iStim)} '.mat'],eegname);eval(['x=' eegname ';' ]);
        [Rxy Rpool] = generate_all_cov(x);
        % load disattend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{disAttend(iStim)} '.mat'],eegname);eval(['y=' eegname ';' ]);
        
        Nsubj = size(y,3);
        % load covariance of attend (normative) group
        [isc_free isc_count] = compute_individual_subject_ISC(x,y,Rpool,Rxy,nComp,gamma,dim,Nsubj);
        iscs{iStim}.free = isc_free;
        iscs{iStim}.count = isc_count;
        clear isc_free isc_count
    end
    save([saveDir],'iscs')
end
