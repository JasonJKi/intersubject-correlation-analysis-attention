clear all
File = [3];
run /home/jason/Experiments/Backwardcount/codes/bdf_reader/eegdetail.m
experiment = {    'freeview'    'fixated'    'audio'    };
processedLoadPath = 'rPCA_processed_eeg/';
% '/home/jason/Experiments/Backwardcount/data/processed/manual/'
savedir='individual_iscs_robust_pca_by_electrode';
mkdir(savedir)

nComp = 8;
nGamma = 9;
grassman = 0;
dim = 9;


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

iFigure = 1
for iFile = [1];
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogchannel=Eogchannel{iFile};
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    subjrmv = Indsubjrmv{iFile};
    
    nStim = length(stims)/2;
    for iStim = 1:nStim;
        load([processedLoadPath experiment{iFile} '_' stimsave{attend(iStim)}],'EEG');eval(['X' '=EEG;' ]);
        % load disattend viewing group
        load([processedLoadPath experiment{iFile} '_' stimsave{disAttend(iStim)}],'EEG');eval(['Y' '=EEG;' ]);
        for iElectrodeSet = 1:nElectrodeSets
            x= X(:,electrodeSets{iElectrodeSet},:);
            y= Y(:,electrodeSets{iElectrodeSet},:);
            [Rxy Rpool] = generate_all_cov(x,0);
            % load covariance of attend (normative) group
            for iGamma =1:nGamma
                [isc_free(:,:,iGamma) isc_count(:,:,iGamma)] = computeIndividualISC(x,y,Rpool,Rxy,nComp,(iGamma*.1),dim);
            end
            iscs{:,iStim}.free{iElectrodeSet} = isc_free;
            iscs{:,iStim}.count{iElectrodeSet} = isc_count;
            clear isc_free isc_count
        end
        save([savedir '/' experiment{iFile}],'iscs')
    end
    clear iscs
end