clear all
eegclean5labpath = '/home/jacek0/Jacek/clean opt/functions/';
addpath([eegclean5labpath 'adminfunc'],...
    [eegclean5labpath 'miscfunc'],...
    [eegclean5labpath 'popfunc'],...
    [eegclean5labpath 'resources'],...
    [eegclean5labpath 'sigprocfunc'],...
    [eegclean5labpath 'studyfunc'],...
    [eegclean5labpath 'timefreqfunc']);

File = [2];
run /home/jason/Experiments/Backwardcount/codes/bdf_reader/eegdetail.m
experiment = {    'freeview'    'fixated'    'audio'    };
loadPath = []; %'/home/jason/Experiments/Backwardcount/data/processed/newmanualmethod/';
nComp = 10;
nGamma = 10;
grassman = 0;
dim = 64;
nCompPlot = 6;
Aall=[]; ISCall=[];
for iFile = File

    
    %File specifications
    attend = Attend{iFile}; disAttend = Disattend{iFile};
    stims=Stims{iFile};stimsave=Stimsave{iFile};
    subjname=Subjname{iFile};eogChannel=Eogchannel{iFile};
    subjrmv = Indsubjrmv{iFile};
    
    nStim = length(stims)/2
    for iStim = 1:nStim
        load([loadPath experiment{iFile} '_' stimsave{attend(iStim)}],'a','isc');
        A(:,:,1) = a; ISC(:,1) = isc;
        % load disattend viewing group
        load([experiment{iFile} '_' stimsave{disAttend(iStim)}],'a','isc');
        A(:,:,2) = a; ISC(:,2) = isc;
        Aall = cat(3,A,Aall); ISCall =cat(2,ISC,ISCall);
    end
       figure(iFile)
       headprojection(Aall,ISCall,nComp)
end
% 
% for iFile = File
%     LoadPath = [processedLoadPath experiment{iFile} '/'];
%     
%     %File specifications
%     attend = Attend{iFile}; disAttend = Disattend{iFile};
%     stims=Stims{iFile};stimsave=Stimsave{iFile};
%     subjname=Subjname{iFile};eogChannel=Eogchannel{iFile};
%     subjrmv = Indsubjrmv{iFile};
%     
%     for iStim = 1:nStim
%         load([experiment{iFile} '_' stimsave{attend(iStim)}],'EEG');
%         [Rxy1(:,:,iStim) Rpool1(:,:,iStim)] =generate_cov(EEG);
%         load([experiment{iFile} '_' stimsave{disAttend(iStim)}],'EEG');
%         [Rxy2(:,:,iStim) Rpool2(:,:,iStim)] =generate_cov(EEG);
%     end
%     [A1 ISC1 w1] = correlated_component(sum(Rxy1,3), sum(Rpool1,3), .5, grassman, dim, nComp,0,1);
%     [A2 ISC2 w2] = correlated_component(sum(Rxy2,3), sum(Rpool2,3), .5, grassman, dim, nComp,0,1);
%     A = concat(3,A1,A2); ISC = concat(3,ISC1,ISC2);
%     figure(iFile)
%     headprojection(A,ISC,nComp)
% end