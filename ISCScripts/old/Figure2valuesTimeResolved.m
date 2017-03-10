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
        
        [T,nChannels,nSubj] = size(x);
        iWind = 0; maxWind = 0
        fs = 256;
        windSize =256*3
        seconds = T/fs
        t = 1;
        while iWind < seconds*t-1*3
            windShift = fs*(iWind/t);
            iWind = iWind +1
            wind = windShift + [1:windSize];
            X = x(wind,:,:);
            Y = y(wind,:,:);
            nSubj=size(y,3);
            % load covariance of attend (normative) group
            [isc_free isc_count] = compute_individual_subject_ISC(X,Y,Rpool,Rxy,nComp,gamma,dim,nSubj);
            iscs(:,:,iWind) = [isc_free isc_count];
        end
        ISCS{iStim,iFile} = iscs;
        clear isc
    end
iFile = iFile + 1
end   

for j = 1:nComp
ISC = sum(isc(nComp,:),1)
for i = 1:size(isc,3)
    Az(i,j) = rocarea(ISC,[ones(1,length(ISCFree)) zeros(1,length(ISCCount))]);
    
end
end