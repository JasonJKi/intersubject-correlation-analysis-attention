%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get All Covariance.
clear all;
homeDir = '/home/jason/AttentionModulation/';
run([homeDir 'eegInfo.m'])
nModule = [1 2 3];

processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
saveDir = 'ISCScripts/ISCValues/'
saveDir = [homeDir saveDir];
covFile =  [processType(1:end-1) '_allCovbyElectrode.mat'];


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
nElec = 8;

i=1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule};
    
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [homeDir eeginfo.eegProcDir processType moduletype '/' stim '.mat'];
        load(eeginDir)
        eeg(:,:,module.subjRmv{iModule})=[];
        for iElec = 1:nElec
            [Rxy_all Rpool_all Rxy Rpool] = generate_cov(eeg(:,electrodeSets{iElec},:),1);
            R{iElec}.xy_all{i} = Rxy_all;
            R{iElec}.pool_all{i} = Rpool_all;
            R{iElec}.xy{i} = Rxy;
            R{iElec}.pool{i} = Rpool;
        end
        i = i+1;
    end
    
end

save([saveDir covFile], 'R')


%% Compute Correlated Components
%load covariance
if ~exist('R')
load([saveDir covFile], 'R');
end 

%correlated component out file
ccFile =  [processType(1:end-1) '_allComponentsbyElectrode.mat'];

%Indx stim groups for generating correlated components and the Forward
%Model
groupIndx{1} = [1 2 3 4 5 6 7 8]
groupIndx{2} = [9 10 11 12 13 14]
nGroup = 2;

%correlated component constants
D = 9; gamma =.5; whitening = 0
for iElec =1:nElec
for iGroup = 1:nGroup;
    Rxy = zeros(D); Rpool=zeros(D);
    for iStim = groupIndx{iGroup}
        Rxy = Rxy + R{iElec}.xy{iStim};
        Rpool = Rpool + R{iElec}.pool{iStim};
    end
    [W{iElec}{iGroup} A{iElec}{iGroup}] = correlated_components(Rxy, Rpool, gamma, D, whitening);
end
end
save([saveDir ccFile], 'W', 'A')

%% compute the ISC.
if ~exist('W')
load([saveDir ccFile],'W')
end

ISCFile =  [processType(1:end-1) '_allISCbyElectrode.mat'];

compIndx = [1 1 2]; nComp = 10;

i = 1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule}
    
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [homeDir eeginfo.eegProcDir processType moduletype '/' stim '.mat'];
        load(eeginDir)
        eeg(:,:,module.subjRmv{iModule}) = [];
        nSubj= size(eeg,3);
        for iElec = 1:nElec
        for iSubj = 1:nSubj
            isc{iElec}(:,iSubj) = concat_matrix_ISC(eeg(:,electrodeSets{iElec},:), ...
            eeg(:,electrodeSets{iElec},:), W{iElec}{compIndx(iModule)},iSubj,nComp);
            disp(iSubj)
        end
        end
        ISC{i} = isc;
        i= i+1;
        clear isc
    end
end
save([saveDir ISCFile], 'ISC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute ISC by leaveoneout.

if ~exist('W')
load([saveDir ccFile],'W')
end

ISCFile =  [processType(1:end-1) '_allISC.mat'];

compIndx = [1 1 2]; nComp = 10;

i = 1
attentionIndx = [ 1 2 5 6 9 10 11;
                    3 4 7 8 12 13 14];

for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule}
    
    for iStim = 1:nStim;
        AttendStim = module.eegStimName{iModule}{attentionIndx(1,iStim)}
        DisattendStim = module.eegStimName{iModule}{attentionIndx(1,iStim)}
        eeginDirAttend = [homeDir eeginfo.eegProcDir processType moduletype '/' stim '.mat'];
        load(eeginDir)
        eeg(:,:,module.subjRmv{iModule}) = [];
        nSubj= size(eeg,3);
        for iSubj = 1:nSubj
        isc(:,iSubj) = concat_matrix_ISC(eeg,eeg,W{compIndx(iModule)},iSubj,nComp);
        disp(iSubj)
        end
        ISC{i} = isc;
        i= i+1;
        clear isc
    end
end
save([saveDir ISCFile], 'ISC')