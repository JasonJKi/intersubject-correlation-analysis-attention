clear all

file = 1:3;
run ../Codes/BDFRead/EEGDetail.m
experiment = { 'Freeview'    'Fixation'    'Audio' };
EEGLoadPath = '../EEG/Processed/';
Version = 'RPCAProcessed_v2/'
EEGLoadPath = [EEGLoadPath Version]
nComp = 20;nGamma = 9;grassman = 0;dim = 64;
saveDir ='individual_iscs_v0'; mkdir(saveDir)

for iFile = 1:3;
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogchannel=Eogchannel{iFile};
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    subjrmv = Indsubjrmv{iFile};
        
    nStim = length(stims)/2;
    for iStim = 1:nStim;
        % load attend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{attend(iStim)} '.mat'],'EEG');eval(['x' '=EEG;' ]);
        [Rxy Rpool] = generate_all_cov(x);
        % load disattend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{disAttend(iStim)} '.mat'],'EEG');eval(['y' '=EEG;' ]);

        % load covariance of attend (normative) group
        for iGamma =1:nGamma
        [isc_free(:,:,iGamma) isc_count(:,:,iGamma)] = compute_individual_subject_ISC(x,y,Rpool,Rxy,nComp,(iGamma*.1),dim);
        end
        iscs{iStim}.free = isc_free;
        iscs{iStim}.count = isc_count;
        clear isc_free isc_count
    end
    save([saveDir '/' experiment{iFile}],'iscs')
end
