% Generate figure 2 optional.

file = 1;
run ../../Codes/BDFRead/EEGDetail.m
experiment = { 'Freeview'    'Fixation'    'Audio' };
EEGLoadPath = '../../EEG/Processed/';
Version = 'RPCAProcessed_v2/'
EEGLoadPath = [EEGLoadPath Version]
nComp = 10;nGamma = 1;grassman = 0;dim = 64; gamma = .5;
 saveDir ='Figure_2_iscs_eogs'; %mkdir(saveDir)

for iFile = file
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogchannel=Eogchannel{iFile};
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    subjrmv = Indsubjrmv{iFile};
        
    nStim = length(stims)/2;
    for iStim = 1:nStim;
        
        % load attend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{attend(iStim)} '.mat'],'eog');eval(['x' '=eog;' ]);
        [Rxy Rpool] = generate_all_cov(x);
        % load disattend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{disAttend(iStim)} '.mat'],'eog');eval(['y' '=eog;' ]);
        
        x(:,:,Indsubjrmv{iFile}) = [];
        y(:,:,Indsubjrmv{iFile}) = [];
        nComp = size(x,2);
        dim = nComp;
        nSubj = size(y,3)
        % load covariance of attend (normative) group
        [isc_free isc_count] = compute_individual_subject_ISC(x,y,Rpool,Rxy,nComp,gamma,dim,nSubj);
        iscs{iFile,iStim}.free = isc_free;
        iscs{iFile,iStim}.count = isc_count;
        clear isc_free isc_count
    end   
end
save([saveDir],'iscs')