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
covFile =  [processType(1:end-1) '_allCov.mat'];

% i=1
% for iModule = nModule
%     nStim = length(module.eegStimName{iModule});
%     eogIndx = module.eogChannel{iModule};
%     moduletype = eeginfo.moduleName{iModule};
%     
%     for iStim = 1:nStim;
%         stim = module.eegStimName{iModule}{iStim}
%         eeginDir = [homeDir eeginfo.eegProcDir processType moduletype '/' stim '.mat'];
%         load(eeginDir)
%         eeg(:,:,module.subjRmv{iModule})=[];
%         [Rxy_all Rpool_all Rxy Rpool] = generate_cov(eeg,1);
%         R.xy_all{i} = Rxy_all;
%         R.pool_all{i} = Rpool_all;
%         R.xy{i} = Rxy;
%         R.pool{i} = Rpool;
%         i = i+1;
%     end
% 
% end
% 
% save([saveDir covFile], 'R')
% 

%% Compute Correlated Components
%load covariance
if ~exist('R')
load([saveDir covFile], 'R');
end 

%correlated component out file
ccFile =  [processType(1:end-1) '_allComponents2.mat'];

%Indx stim groups for generating correlated components and the Forward
%Model
% groupIndx{1} = [1 2 3 4 5 6 7 8]
% groupIndx{2} = [9 10 11 12 13 14]
groupIndx{1} = [1 2 4 5]
groupIndx{2} = [9 10 11]

nGroup = 2;

%correlated component constants
D = 64; gamma =.5; whitening = 0

for iGroup = 1:nGroup;
    Rxy = zeros(64); Rpool=zeros(64);
    for iStim = groupIndx{iGroup}
        Rxy = Rxy + R.xy{iStim};
        Rpool = Rpool + R.pool{iStim};
    end
    [W{iGroup} A{iGroup}] = correlated_components(Rxy, Rpool, gamma, D, whitening);
end

save([saveDir ccFile], 'W', 'A')

%% compute the ISC.
if ~exist('W')
load([saveDir ccFile],'W')
end

ISCFile =  [processType(1:end-1) '_allISC.mat'];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute ISC by leaveoneout.

if ~exist('W')
load([saveDir ccFile],'W')
end

ISCFile =  [processType(1:end-1) '_allISCPredictor.mat'];

compIndx = [1 1 2]; nComp = 10;

attentionIndx = {[1 2; 3 4], [1 2; 3 4],[1 2 3; 4 5 6]};
storeIndx = [1 2 5 6 9 10 11;
             3 4 7 8 12 13 14];
i = 1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule};
    
    for iStim = 1:nStim/2;
        AttendStim = module.eegStimName{iModule}{attentionIndx{iModule}(1,iStim)}
        DisattendStim = module.eegStimName{iModule}{attentionIndx{iModule}(2,iStim)}
        eeginDirAttend = [homeDir eeginfo.eegProcDir processType moduletype '/' AttendStim '.mat'];
        load(eeginDirAttend)
        eeg(:,:,module.subjRmv{iModule}) = [];
        x = eeg;
        eeginDirDisAttend = [homeDir eeginfo.eegProcDir processType moduletype '/' DisattendStim '.mat'];
        load(eeginDirDisAttend)
        eeg(:,:,module.subjRmv{iModule}) = [];
        y = eeg;
        clear eeg
        nSubj= size(y,3);
        parfor iSubj = 1:nSubj
        iscA(:,iSubj) = concat_matrix_ISC(x,x,W{compIndx(iModule)},iSubj,nComp);
        iscD(:,iSubj) = concat_matrix_ISC(y,x,W{compIndx(iModule)},iSubj,nComp);
        disp(iSubj)
        end
        ISC{storeIndx(1,i)} = iscA;
        ISC{storeIndx(2,i)} = iscD;
        i= i+1;
        clear isc 
    end
end
save([saveDir ISCFile], 'ISC')