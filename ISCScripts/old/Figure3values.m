% Generate figure 2.
clear all
file = 1:3;
run ../../Codes/BDFRead/EEGDetail.m
experiment = { 'Freeview'    'Fixation'    'Audio' };
EEGLoadPath = '../../EEG/Processed/';
eegname ='eeg'
Version = 'RPCAProcessed_v4/'
EEGLoadPath = [EEGLoadPath Version]
nComp = 8;nGamma = 1;grassman = 0;dim = 9; gamma = .9;
 saveDir ='Figure_3_iscs_v2'; %mkdir(saveDir)

electrodeSets = {
    [1  2  3  33 34 35 36 37 38];...
    [7  6  5  4  38 39 40 41 42];...
    [8  9  10 11 47 46 45 44 43];...
    [15 14 13 12 48 49 50 51 52];...
    [16 17 18 19 32 56 55 54 53];...
    [23 22 21 20 31 57 58 59 60];...
    [25 26 27 29 30 31 62 63 64];...
    [24 25 27 28 29 30 61 62 64];...
    };

nElectrodeSets = size(electrodeSets,1);

for iFile = file;
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogchannel=Eogchannel{iFile};
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    subjrmv = Indsubjrmv{iFile};
    
    nStim = length(stims)/2;
    for iStim = 1:nStim;
         % load attend viewing group
        % load attend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{attend(iStim)} '.mat'],eegname);eval(['x=' eegname ';' ]);
        [Rxy Rpool] = generate_all_cov(x);
        % load disattend viewing group
        load([EEGLoadPath experiment{iFile} '/' stimsave{disAttend(iStim)} '.mat'],eegname);eval(['y=' eegname ';' ]);
  
        for iElectrodeSet = 1:nElectrodeSets
            X= x(:,electrodeSets{iElectrodeSet},:);
            Y= y(:,electrodeSets{iElectrodeSet},:);
            nSubj = size(Y,3);
            [Rxy Rpool] = generate_all_cov(X);
            % load covariance of attend (normative) group
            [isc_free isc_count] = compute_individual_subject_ISC(X,Y,Rpool,Rxy,nComp,gamma,dim,nSubj);
            iscs{iFile,iStim}.free{iElectrodeSet} = isc_free;
            iscs{iFile,iStim}.count{iElectrodeSet} = isc_count;
        end
    end
end
save(saveDir,'iscs')


map = [0, 0, 0.3
    0, 0, 0.4
    0, 0, 0.5
    0, 0, 0.6
    0, 0, 0.8
    0, 0, 1.0];

figure
k =surf(peaks)
colormap(map)
