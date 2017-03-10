%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get All Covariance.
clear all;
homeDir = '/home/jason/AttentionModulation/';
run([homeDir 'eegInfo.m'])
nModule = [1 2 3];

processIndx = 1 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
saveDir = 'ISCScripts/ISCValues/'
saveDir = [homeDir saveDir];
covFile =  [processType(1:end-1) '_allCovEOG.mat'];

i=1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule};
    
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eoginDir = [homeDir eeginfo.eegProcDir processType moduletype '/' stim '.mat'];
        load(eoginDir,'eog')
        eog(:,:,module.subjRmv{iModule})=[];
        for iElec = 1
            [Rxy_all Rpool_all Rxy Rpool] = generate_cov(eog,1);
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
ccFile =  [processType(1:end-1) '_allComponentsEOG.mat'];

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

ISCFile =  [processType(1:end-1) '_allISCEOG.mat'];

compIndx = [1 1 2]; nComp = 9;

i = 1
for iModule = nModule
    nStim = length(module.eogStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eoginfo.moduleName{iModule}
    
    for iStim = 1:nStim;
        stim = module.eogStimName{iModule}{iStim}
        eoginDir = [homeDir eoginfo.eogProcDir processType moduletype '/' stim '.mat'];
        load(eoginDir)
        eog(:,:,module.subjRmv{iModule}) = [];
        nSubj= size(eog,3);
        parfor iElec = 1:nElec
            for iSubj = 1:nSubj
            isc{iElec}(:,iSubj) = concat_matrix_ISC(eog(:,electrodeSets{iElec},:), ...
            eog(:,electrodeSets{iElec},:), W{iElec}{compIndx(iModule)},iSubj,nComp);
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

ISCFile =  [processType(1:end-1) '_allISCbyEOGPredictor.mat'];

compIndx = [1 1 2]; nComp = 9;


attentionIndx = {[1 2; 3 4], [1 2; 3 4],[1 2 3; 4 5 6]};
storeIndx = [1 2 5 6 9 10 11;
             3 4 7 8 12 13 14];
i = 1

for iModule = nModule
    nStim = length(module.eogStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eoginfo.moduleName{iModule}
    
    for iStim = 1:nStim/2;
        AttendStim = module.eogStimName{iModule}{attentionIndx{iModule}(1,iStim)}
        DisattendStim = module.eogStimName{iModule}{attentionIndx{iModule}(2,iStim)}
        eoginDirAttend = [homeDir eoginfo.eogProcDir processType moduletype '/' AttendStim '.mat'];
        eoginDirDisttend = [homeDir eoginfo.eogProcDir processType moduletype '/' DisattendStim '.mat'];
        load(eoginDirAttend)
        eog(:,:,module.subjRmv{iModule}) = [];
        x = eog;
        load(eoginDirDisttend)
        eog(:,:,module.subjRmv{iModule}) = [];
        y = eog;
        clear eog
        nSubj= size(y,3);
        parfor iElec = 1:nElec
            for iSubj = 1:nSubj
            iscA{iElec}(:,iSubj) = concat_matrix_ISC(x(:,electrodeSets{iElec},:), ...
            x(:,electrodeSets{iElec},:), W{iElec}{compIndx(iModule)},iSubj,nComp);
            iscD{iElec}(:,iSubj) = concat_matrix_ISC(y(:,electrodeSets{iElec},:), ...
            x(:,electrodeSets{iElec},:), W{iElec}{compIndx(iModule)},iSubj,nComp);
        disp(iSubj)
        end
        end
        ISC{storeIndx(1,i)} = iscA;
        ISC{storeIndx(2,i)} = iscD;
        i= i+1;
        clear isc
    end
end
save([saveDir ISCFile], 'ISC')