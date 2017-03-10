% Get All Covariance.
clear all;
homeDir = '/home/jason/Repository/AttentionModulation/';
nModule = [1 2];

saveDir = 'Data/EEG/Processed/OrderReversed_v2/'
saveDir = [homeDir saveDir];
stimNames={'bng_1' 'bng_2' 'gbu_1' 'gbu_2' 'pieOrg_1' 'pieOrg_2'}
nStim=length(stimNames);

i=1
for iStim = 1:nStim;
    eeginDir = [saveDir stimNames{iStim} '.mat'];
    load(eeginDir)
    [Rxy Rpool] = generate_cov(eeg,1);
    R.xy{i} = Rxy;
    R.pool{i} = Rpool;
    i = i+1
end

groupIndx{1} = [1 2 3 4];
groupIndx{2} = [5 6];
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

% check the components
figure(1);headprojection(A{1},3,'hor',0);
figure(2);headprojection(A{2},3,'hor',0);

nComp=10;
compIndx=[1 1];
filter.status=0;

i=1; 
for iModule=1:2
       
    for iStim = groupIndx{iModule};
        eeginDir = [saveDir stimNames{iStim} '.mat'];
        load(eeginDir)
        nSubj= size(eeg,3);
        for iSubj = 1:nSubj
            isc(:,iSubj) = concat_matrix_ISC(eeg,eeg,W{compIndx(iModule)},iSubj,nComp);
            disp(iSubj)
        end
        ISC{i} = isc;
        clear isc
            i= i+1;
    end
end

save('ISC_reverseOrder2.mat','ISC')

ISC2= ISC;
load('RPCA_allISC.mat')
ISC1= ISC;
ISC=[ISC1 ISC2];
save('ISC_allConditions2.mat','ISC')
load(['ISC_allConditions2.mat'],'ISC')


i = 1
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule}
    
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [homeDir eeginfo.eegProcDir processtype moduletype '/' stim '.mat'];
        load(eeginDir);
        
        eeg(:,:,module.subjRmv{iModule}) = [];
        if i == 1
            eeg(:,:,end) = [];
        end
        nSubj= size(eeg,3);
        for iSubj = 1:nSubj
            power{i}(:,iSubj)= mean(abs(filter(b,a,eeg(:,:,iSubj)*W{compIndx(iModule)}(:,1:nComp)).^2),1);
             power{i}(:,iSubj) =  power{i}(:,iSubj)./mean(abs(eeg(:,:,iSubj)*W{compIndx(iModule)}(:,1:nComp)).^2,1)';
            disp(iSubj)
        end
        i= i+1;
    end
end

save([homeDir currentDir 'alpha_power_gabor2_theta.mat'], 'power')

%% figure 6a- ATTENTION COMPARISON (BYD,GBU,PM,PMsc,Jpn) alpha power


%create alpha band filter
n = 7
fsRef =256
[b,a] = butter(n,[8 12]/fsRef*2,'bandpass');
freqz(b,a,fsRef,fsRef)
nComp = 3

% check the components
headprojection(A{1},3,'hor',0);
figure(2)
headprojection(A{2},3,'hor',0);

fs = fsRef; % sampling rate in Hz
fc = 5; % center frequency in Hz
Q = 1;  % Q-factor = f/df;
df = 5; %fc/Q; % bandwidth in Hz - wider (used to be 3)
dt = 1/df;
t = (-3*dt*fs:3*dt*fs)'/fs;
b = 1/sqrt(pi/2)/fs/dt*exp(-t.^2/2/dt^2).*exp(sqrt(-1)*2*pi*fc*t);
a = 1;
freqz(b,a,fsRef,fsRef)

i = 1
compIndx=[1 2];
for iModule=1:2
for iStim = 1:nStim;
    eeginDir = [saveDir stimNames{iStim} '.mat'];
    load(eeginDir)
    nSubj= size(eeg,3);
    nSamples=size(eeg,1);
    nWindows=mod(nSamples,fsRef);
    for iSubj = 1:nSubj
        for iW=1:nWindows
            window=(1-iW)*fsRef+1:(iW)*fsRef;
            powerWindow(iW)=mean(abs(filter(b,a,eeg(window,:,iSubj)*W{compIndx(iModule)}(:,1:nComp)).^2),1);
        end
        disp(iSubj)
        power_{i}(:,iSubj)=powerWindow
    end
    i= i+1;
end
end


% Stim labels, and indexes for plotting
modality = [1 1 1 1 2 2 ];
attention = [1 2 1 2 1 2 1 2 1 2 1 2 1 2];
stimLabel = [1 1 2 2 3 3];
visualIndx = [1 1 1 1 0 0];
xAxis = [1 2 3 4 5 6];
narrative = [1 1 1 1 1 1];

% initialize variables
powerAll = []; powerSumAll = []; powerByCompAll =[]; powerRelativeContributionAll = [];
powerRelativeContributionMeanAll=[]; 
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
powerCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];

% set components of interest
compInt = 1:3 % components of interest
nComp = length(compInt) % measure first Ncomps

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6];[.5 .2 .8]};
stimIndx = [1 1 1 1 1 1];
subjgroup = [1 1 1 1 1 1];
stimOrder = [1 2 3 4 5 6];
nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;
f =figure(1);
clf
subplot(2,9,[1:7])
for iStim = 1:nStim
    power_ = power{stimOrder(iStim)}(compInt,subjInt);
    
    powerSum = sum(power_,1);
    powerSumMean(iStim) = mean(sum(power_,1));
    
    if iStim == 1; powerSum(end) =[];end
    %labels
    nSubj=length(powerSum)
    subj = [1:nSubj]+((subjgroup(iStim))*size(powerSum,2)) ;
    attentionIndx = repmat(attention(iStim), nSubj,1);
    modalityIndx = repmat(modality(iStim), nSubj,1);
    narrativeIndx = repmat(narrative(iStim),nSubj,1);
    stimIndx = repmat(stimLabel(iStim),nSubj,1);
    
    if stimIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        powerSumAll = [powerSumAll; powerSum'];
        narrativeIndxAll = [narrativeIndxAll; narrativeIndx];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx];
        subjLabelAll = [subjLabelAll; subj'];
        stimIndxAll = [stimIndxAll; stimIndx];
        % plottings
        hold on;
        if attention(iPlot) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),mean(powerSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        %                 plot(xAxis(iPlot),powerSum,'.');
        iPlot = iPlot +1
    end
end

videos = {'BYD' 'GBU' 'PM'};
mode = {'audiovisual(free)' 'audiovisual(constrained)' 'audio only'} 
set(gca, 'Xtick',  [2.5 8 13.5 ], 'Xticklabel', mode);
% plot properties
p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
l1=legend([h([1 3 5]) p1 p2],[videos 'attend', 'count'])
rect = [0.76, 0.75, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
axis([min(xAxis)-1 max(xAxis)+1 .0 .3])
ylabel('alpha power')
box Off
% saveas(f,'figure6a_attendVsDisattend','epsc')

%% compute significance figure 1
stimIndx = [1 3 5; 2 4 6];
nStim = size(stimIndx,2);
for iStim = 1:3
    attend = sum(power{stimIndx(1,iStim)}(compInt,subjInt));
    disattend = sum(power{stimIndx(2,iStim)}(compInt,subjInt));
    power_ = [attend disattend]';
    attentionIndx= [ones(length(attend),1) ;zeros(length(disattend),1)];
    [h1(iStim) p(iStim)] = ttest2(attend,disattend);
%     [p_{iStim} t]= anovan(power_,attentionIndx,'display','off');
%     F(iStim) = t{2,6};
end

nesting = [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 1 0 0];
model = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 1 0; 0 0 0 1 1];
[P1,T1,STATS1,TERMS1] = anovan(powerSumAll,[subjLabelAll modalityIndxAll attentionIndxAll stimIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1 4], ...
    'varname', {'subject' 'modality' 'attention' 'stim'}, 'display', 'on')

%{
Notes:
UPDATED April 3rd
figure 1 - add n.s. with text height issue corrected 
figure 2 - add component significance line to each bar.
         - add n.s. with text height issue corrected
figure 3 - change significance line
figure 4 - add only visual free byd, gbu, PM
figure 5 - add n.s. with text height issue corrected
%}

% Generate figures and Stats of All ISC
homeDir = '/home/jason/Repository/AttentionModulation/';
run([homeDir 'Codes/eegInfo.m'])
nModule = [1 2 3];

processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
iscFile =  [processType(1:end-1) '_allISC.mat'];
saveDir = 'Data/ISCValues/'
saveDir = [homeDir saveDir];
load([saveDir iscFile]);

stims = {'BYD attend'
    'GBU attend'
    'BYD Count'
    'GBU Count'
    'BYD constraint attend'
    'GBU constraint attend'
    'BYD constraint count'
    'GBU constraint count'
    'PM attend'
    'PM scrm attend'
    'Jpn attend'
    'PM count'
    'PM scrm count'
    'Jpn Count'};

videos = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};

% ISC1= ISC;
% load('ISC_reverseOrder.mat')
% ISC2=ISC;
% ISC=[ISC1 ISC2];
% save('ISC_allConditions.mat','ISC')
load([saveDir 'ISC_allConditions.mat'],'ISC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 7 - ATTENTION COMPARISON (BYD,GBU,PM,PMsc,Jpn)
% Stim labels, and indexes for plotting
modality = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1 1 2 2];
stimLabel = [1 1 2 2 1 1 2 2 3 3 4 4 5 5 1 1 2 2 3 3];
visualIndx =[1 1 1 1 2 2 2 2 0 0 0 0 0 0 1 1 1 1 0 0];
viewingOrder=  [1 1 1 1 1 2 1 2 1 1 1 2 1 2 1 2 1 2 1 2];
xAxis = [1:6 8:13];

% initialize variables
iscAll = []; iscSumAll = []; iscByCompAll =[]; iscRelativeContributionAll = [];
iscRelativeContributionMeanAll=[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];orderIndxAll=[];
iscCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];

% set components of interest
compInt = 1:3 % components of interest
nComp = length(compInt) % measure first Ncomps 

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; };
figure(1); clf; sp1= subplot(5,5,[1 2 3 6 7 8 11 12 13]); 

stimIndx = [1 1 1 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0];
subjgroup = [1 1 1 1 0 0 0 0 2 2 2 2 2 2 3 3 3 3 3 3 ];
stimOrder = [1 16 2 18 5 7 6 8 9 20 10 13 11 14 16 15 18 17 20 19];
attention = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 2 1 2 1 2 1];

nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;
for iStim = 1:nStim
    isc_ = ISC{stimOrder(iStim)}(compInt,:);
    isc__ = ISC{stimOrder(iStim)}(1:10,:);
    
    iscRelativeContribution = isc_./repmat(sum(abs(isc_),1),length(compInt),1);    
    iscRelativeContributionAll = [iscRelativeContributionAll iscRelativeContribution];
        
    iscRelativeContributionMean = mean(iscRelativeContribution,2);
    iscRelativeContributionMeanAll(:,iStim) = iscRelativeContributionMean;
        
    iscSum = sum(isc_,1);
    iscByComp(:,iStim) = mean(isc__,2);
    iscSumMean(iStim) = mean(sum(isc_,1));
   
    if iStim == 1; iscSum(end) =[];end
    %labels
    nSubj=length(iscSum)
    subj = [1:nSubj]+((subjgroup(iStim))*size(iscSum,2)) ;
    attentionIndx = repmat(attention(iStim), nSubj,1);
    modalityIndx = repmat(modality(iStim), nSubj,1); 
    orderIndx =repmat(viewingOrder(iStim), nSubj,1); 
    if stimIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        iscSumAll = [iscSumAll; iscSum'];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx]; 
        subjIndxAll = [subjIndxAll; subj'];
        orderIndxAll = [orderIndxAll; orderIndx]
        % plottings
        hold on;
        if attention(iStim) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),mean(iscSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        plot(xAxis(iPlot),iscSum,'.');
        n{iPlot} = num2str(length(iscSum));
        iPlot = iPlot +1
    end
end

p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
l1=legend([p1 p2],'attend', 'count')
rect = [0.475, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
delete([p1 p2]);
ylabel('ISC') ;
plotIndx = [3.5 10.5];
set(gca, 'Xtick', plotIndx,'xticklabel',{'ordering independent' 'count-attend(reversed viewing order)'})
t1=title('a');set(t1,'position',[-.05 1])
saveas(figure(1),'figure7b','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%
nesting = [0 1 1 ; 0 0 0 ; 0 0 0]
labelIndx=[subjIndxAll modalityIndxAll attentionIndxAll];
labels={'subject' 'modality' 'attention' }

nesting = [];nesting = [0 1 1 ; 0 0 0 ; 0 0 0];
labelIndx=[subjIndxAll attentionIndxAll orderIndxAll];
labels={ 'subject' 'attention' 'order'}

[P1,T1,STATS1,TERMS1] = anovan(iscSumAll,labelIndx, ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',1, ...
    'varname', labels, 'display', 'on')

%% figure 7 - ATTENTION COMPARISON (BYD,GBU,PM,PMsc,Jpn)
% Stim labels, and indexes for plotting
modality =      [1 1 1 1 2 2 1 1 1 1 2 2];
stimLabel =     [1 1 2 2 3 3 1 1 2 2 3 3];
visualIndx =    [1 1 1 1 0 0 1 1 1 1 0 0];
viewingOrder=   [1 1 1 1 1 1 2 2 2 2 2 2];
xAxis = [1:6 8:13];

% initialize variables
iscAll = []; iscSumAll = []; iscByCompAll =[]; iscRelativeContributionAll = [];
iscRelativeContributionMeanAll=[];
subjIndxAll= [];stimIndxAll=[];modalityIndxAll= [];orderIndxAll=[];
iscCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];

% set components of interest
compInt = 1:3 % components of interest
nComp = length(compInt) % measure first Ncomps 

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; };
figure(1); clf; sp1= subplot(5,5,[1 2 3 6 7 8 11 12 13]); 

plotIndx  = [1 0 1 0 1 0 1 0 1 0 1 0]; 
subjgroup = [1 2 1 2 3 2 1 2 1 2 2 3];
stimOrder = [1 16 2 18 9 20 15 3 17 4 19 12];
attention = [1 2 1 2 1 2 1 2 1 2 1 2];

nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;

for iStim = 1:nStim
    isc_ = ISC{stimOrder(iStim)}(compInt,:);
    isc__ = ISC{stimOrder(iStim)}(1:10,:);
    
    iscRelativeContribution = isc_./repmat(sum(abs(isc_),1),length(compInt),1);    
    iscRelativeContributionAll = [iscRelativeContributionAll iscRelativeContribution];
        
    iscRelativeContributionMean = mean(iscRelativeContribution,2);
    iscRelativeContributionMeanAll(:,iStim) = iscRelativeContributionMean;
        
    iscSum = sum(isc_,1);
    iscByComp(:,iStim) = mean(isc__,2);
    iscSumMean(iStim) = mean(sum(isc_,1));
   
    if iStim == 1; iscSum(end) =[];end
    %labels
    nSubj=length(iscSum);
    subj = [1:nSubj]+((subjgroup(iStim))*size(iscSum,2)) ;
    attentionIndx=repmat(attention(iStim), nSubj,1);
    modalityIndx=repmat(modality(iStim), nSubj,1); 
    orderIndx=repmat(viewingOrder(iStim), nSubj,1); 
    stimIndx=repmat(stimLabel(iStim),nSubj,1);
    if plotIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        iscSumAll = [iscSumAll; iscSum'];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx]; 
        subjIndxAll = [subjIndxAll; subj'];
        orderIndxAll = [orderIndxAll; orderIndx]
        stimIndxAll = [stimIndxAll; stimIndx]
% plottings
        hold on;
        if attention(iStim) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),mean(iscSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        plot(xAxis(iPlot),iscSum,'.');
        n{iPlot} = num2str(length(iscSum));
        iPlot = iPlot +1
    end
end


p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
l1=legend([h([1 3 5]) p1 p2 ],{'BYD','GBU','PM','attend','count'})
rect = [0.5, 0.70, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
delete([p1 p2]);
ylabel('ISC') ;
plotIndx = [3.5 10.5];
set(gca, 'Xtick', plotIndx,'xticklabel',{'1st viewing' '2nd viewing'})
t1=title('a');set(t1,'position',[-.05 1])

%pair wise comparison

% compute significance 
stimIndx = [1 2 9 15 17 19; 16 18 20 3 4 12];
for iStim = 1:size(stimIndx,2)
    attend = sum(ISC{stimIndx(1,iStim)}(compInt,:));
        if iStim==1;attend(end)=[];end
    disattend = sum(ISC{stimIndx(2,iStim)}(compInt,:));
    isc = [attend disattend]';
    subjIndx=[1:length(attend) 1:length(disattend)]';
    attentionIndx= [ones(length(attend),1) ;zeros(length(disattend),1)];
    [h1(iStim),p(iStim)] = ttest2(attend,disattend);
    [p_ t]= anovan(isc,[attentionIndx],....
             'display','off');
%     P(iStim)=p_;
    F(iStim) = t{2,6};
end

linestyle = {'-k' '-k' '-k' '-k' '-k' '-k'};
yheight = [0 0 0 0 0 0];
pBonf = bonf_holm(p);
pBonf(pBonf > .05) = nan;
pFdr =mafdr(p,'BHFDR',1);
pFdr(pFdr > .05) = nan;
p(p > .05) = nan
sigstarCustom({[1 2],[3 4], [5 6], [8 9], [10 11],[12 13]},pBonf,[],linestyle,yheight)

% nesting = [0 1 1 ; 0 0 0 ; 0 0 0];
nesting = [];
labelIndx=[orderIndxAll stimIndxAll];
labels={'order' 'stimulus'}

[P1,T1,STATS1,TERMS1] = anovan(iscSumAll,labelIndx, ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[], ...
    'varname', labels, 'display', 'on')
% saveas(figure(3),'stats','epsc')

stimIndx =[1 2 9 16 18 20;15 17 19 3 4 12];
for iStim = 1:size(stimIndx,2)
    attend = sum(ISC{stimIndx(1,iStim)}(compInt,:));
        if iStim==1;attend(end)=[];end
    disattend = sum(ISC{stimIndx(2,iStim)}(compInt,:));
    isc = [attend disattend]';
    subjIndx=[1:length(attend) 1:length(disattend)]';
    attentionIndx= [ones(length(attend),1) ;zeros(length(disattend),1)];
    [h1(iStim),p(iStim)] = ttest2(attend,disattend);
    [p_ t]= anovan(isc,[attentionIndx],....
     'display','off');
    P(iStim)=p_;
    F(iStim) = t{2,6};
end

linestyle = {'-g' '-g' '-g' '-b' '-b' '-b'};
yheight = [.163 .17 .177 .19 .202 .214];
pBonf = bonf_holm(P);
pBonf(pBonf > .05) = nan;
pFdr =mafdr(p,'BHFDR',1);
pFdr(pFdr > .05) = nan;
p(p > .05) = nan
sigstarCustom({[1 8],[3 10], [5 12], [2 9], [4 11],[6 13]},pBonf,[],linestyle,yheight)

axis([0 14 0 .215])

saveas(figure(1),['../../Data/' 'figure6'],'epsc')


% 
% figure(2);clf
% subplot(2,3,[1 2])
% plotIndx = [1 1 2 2 3 3 1 1 2 2 3 3 1 1 2 2 3 3 ];
% indxViewing=[1 3 5; 2 4 6];
% barColIndx = [1 2 3]
% barCol = [.3 .9 .9; 0 .6 .6;.5 .2 .8];
% iPlot = 1;
% i = 1
% 
% subjInt=[1:2 4 6:20];;
% 
% for iComp = 1:3
%     for iStim = 1:3
%         attendISC = ISC2{indxViewing(1,iStim)}(iComp,subjInt)';
%         disattendISC = ISC2{indxViewing(2,iStim)}(iComp,subjInt)';
%         attendISCMean = mean(attendISC,1);
%         disattendISCMean = mean(disattendISC,1);
% %         attendISCChanceMean = iscChanceMean(iComp,indxViewing(1,iStim)); 
% %         disattendISCChanceMean = iscChanceMean(iComp,indxViewing(2,iStim));
%          
%         [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendISC,disattendISC); 
%         % bar graph
%         hold on
%         if (iComp == 3 & iStim == 3)
% %             k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%             h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
%             iPlot = iPlot +1
% %             k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%             h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
%         else
%             h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
% %             k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%             iPlot = iPlot +1
%             h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
% %             k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%         end
%     end
% end
% 
%         
% figure(2);clf
% subplot(2,3,[1 2])
% plotIndx = [1 1 2 2 3 3 1 1 2 2 3 3 1 1 2 2 3 3 ];
% indxViewing=[1 3 5; 2 4 6];
% barColIndx = [1 2 3]
% barCol = [.3 .9 .9; 0 .6 .6;.5 .2 .8];
% iPlot = 1;
% i = 1
% 
% subjInt=[1:2 4 6:20];;
% 
% for iComp = 1:3
%     for iStim = 1:3
%         attendISC = ISC2{indxViewing(1,iStim)}(iComp,subjInt)';
%         disattendISC = ISC2{indxViewing(2,iStim)}(iComp,subjInt)';
%         attendISCMean = mean(attendISC,1);
%         disattendISCMean = mean(disattendISC,1);
% %         attendISCChanceMean = iscChanceMean(iComp,indxViewing(1,iStim)); 
% %         disattendISCChanceMean = iscChanceMean(iComp,indxViewing(2,iStim));
%          
%         [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendISC,disattendISC); 
%         % bar graph
%         hold on
%         if (iComp == 3 & iStim == 3)
% %             k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%             h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
%             iPlot = iPlot +1
% %             k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%             h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
%         else
%             h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
%                 'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
% %             k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
% %                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
%             iP
%         %graph property
%                 
%         yHeightMax = max([attendISCMean disattendISCMean]);
% %         yHeightMax = max([attendISCMean disattendISCMean attendISCChanceMean disattendISCChanceMean]);
%         yHeight(i) = yHeightMax;
%         iPlot = iPlot +1
%         i = i+1;
%     end
% end
% 
% p =mystack(cell2mat(pValue'))
% pBonf = bonf_holm(p);
% pBonf(pBonf > .05) = nan
% pFdr =mafdr(p,'BHFDR',2);
% pFdr(pFdr > .05) = nan;
% p(p > .05) = nan
% lineStyle = {'-k', '-k', '-k', '-k', '-k', '-k','-k', '-k', '-k'};
% sigstarCustom({[1 2], [3 4], [5 6], [8 9], [10 11], [12 13], [15 16], [17 18],[19 20]},pBonf,[],lineStyle,yHeight+.001);
% set(gca, 'Xtick', [3.5 10.5 17.5], 'Xticklabel', {'C1', 'C2' ,'C3'});
% p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
% p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
% axis([0 21 0 .035])
% l1=legend([h([2 4 6]) p2 p1], {'BYD' 'GBU' 'PM' 'count' 'attend'})
% delete([p1 p2]);
% ylabel('ISC')
% box Off
% 
% saveas(figure(2),'Figure6','epsc')
% 
% 
% 
% %create alpha band filter
% n = 7
% fsRef =256
% [b,a] = butter(n,[8 12]/fsRef*2,'bandpass');
% freqz(b,a,fsRef,fsRef)
% nComp = 3
% 
% % check the components
% figure(1)
% headprojection(A{1},3,'hor',0);
% figure(2)
% headprojection(A{2},3,'hor',0);
% 
% 
% fs = fsRef; % sampling rate in Hz
% fc = 5; % center frequency in Hz
% Q = 1;  % Q-factor = f/df;
% df = 5; %fc/Q; % bandwidth in Hz - wider (used to be 3)
% dt = 1/df;
% t = (-3*dt*fs:3*dt*fs)'/fs;
% b = 1/sqrt(pi/2)/fs/dt*exp(-t.^2/2/dt^2).*exp(sqrt(-1)*2*pi*fc*t);
% a = 1;
% freqz(b,a,fsRef,fsRef)

% attention =[1 2 9 10 11; 3 4 12 13 14];
% stimOrder= [1 3 6 8 10; 2 4 7 9 11];
% azIndx = [1 2 4 5 6];
% 
% subplot(2,2,2)
% barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
% %Az Significance.
% Nshuffles = 10000;
% for iStim =1:5
%     iscAttend = sum(ISC{attention(1,iStim)}(compInt,:));
%     iscCount = sum(ISC{attention(2,iStim)}(compInt,:));
%     label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
%     az(iStim) = rocarea([iscAttend iscCount],label);
%     lengthlabel = length(label);
%     for iShuffle = 1:Nshuffles
%         az_(iShuffle) = rocarea([iscAttend iscCount], label(randperm(lengthlabel)));
%     end
%     p = length(find(az_ > az(iStim)))/Nshuffles;
%     i = 10
%     rnd = 1
%     while rnd
%         if p == 0
%             ps = num2str(.0001,'%.4f');
%             rnd =0
%         elseif p < 1/10^(i)
%             i
%             p_=floor(p*10^(i+2))/10^(i+2);
%             ps= num2str(p_);
%             rnd = 0;
%         end
%         i= i-1
%     end
%     P{iStim} = num2str(ps,'%.4f');
% end
% 
% cdfAz = cdf(makedist('Normal',mean(az_),std(az_)),az_)
% sigPvalue = min(az_(find(cdfAz> .99))) 
% 
% % shade area of insignificance
% x=[0 (nStim+2)]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
% shadedplot([0 7], y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on
% 
% for iStim =1:5
%     hold on
%     iscAttend = sum(ISC{attention(1,iStim)}(compInt,:),1);
%     iscCount = sum(ISC{attention(2,iStim)}(compInt,:),1);
%     label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
%     az(iStim) = rocarea([iscAttend iscCount], label);
%     bar(azIndx(iStim), az(iStim),0.5,'FaceColor', barcol{iStim}, ...
%          'EdgeColor',barcol{iStim}*.02,'LineWidth',1,'BarWidth',.75)
% end
% 
% plot([0 7],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1,'Color', [.25 .25 .25])
% ylabel('Az');stims = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
% set(gca, 'Xtick', azIndx ,'xticklabel', stims, 'Ytick', [.5 .75 1], 'Yticklabel', {'0.5' '.75' '1'});
% axis([0 7 .5 1])
% saveas(figure(3), 'figure3', 'epsc')
% 
% 
% %Az Significance.
% Nshuffles = 10000;
% for iStim =1:5
%     iscAttend = sum(ISC{attention(1,iStim)}(compInt,:));
%     iscCount = sum(ISC{attention(2,iStim)}(compInt,:));
%     label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
%     az(iStim) = rocarea([iscAttend iscCount],label);
%     lengthlabel = length(label);
%     for iShuffle = 1:Nshuffles
%         az_(iShuffle) = rocarea([iscAttend iscCount], label(randperm(lengthlabel)));
%     end
%     p = length(find(az_ > az(iStim)))/Nshuffles;
%     i = 10
%     rnd = 1
%     while rnd
%         if p == 0
%             ps = num2str(.0001,'%.4f');
%             rnd =0
%         elseif p < 1/10^(i)
%             i
%             p_=floor(p*10^(i+2))/10^(i+2);
%             ps= num2str(p_);
%             rnd = 0;
%         end
%         i= i-1
%     end
%     P{iStim} = num2str(ps,'%.4f');
% end
% 
% cdfAz = cdf(makedist('Normal',mean(az_),std(az_)),az_)
% sigPvalue = min(az_(find(cdfAz> .99))) 
% 
% % shade area of insignificance
% x=[0 (nStim+2)]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
% shadedplot([0 7], y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on
% 
% for iStim =1:5
%     hold on
%     iscAttend = sum(ISC{attention(1,iStim)}(compInt,:),1);
%     iscCount = sum(ISC{attention(2,iStim)}(compInt,:),1);
%     label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
%     az(iStim) = rocarea([iscAttend iscCount], label);
%     bar(azIndx(iStim), az(iStim),0.5,'FaceColor', barcol{iStim}, ...
%          'EdgeColor',barcol{iStim}*.02,'LineWidth',1,'BarWidth',.75)
% end
% 
% plot([0 7],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1,'Color', [.25 .25 .25])
% ylabel('Az');stims = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
% set(gca, 'Xtick', azIndx ,'xticklabel', stims, 'Ytick', [.5 .75 1], 'Yticklabel', {'0.5' '.75' '1'});
% axis([0 7 .5 1])
% saveas(figure(3), 'figure3', 'epsc')
% 
% nesting = [];
% labelIndx=[ attentionIndxAll orderIndxAll];
% labels={'attention' 'order'}
% 
% [P1,T1,STATS1,TERMS1] = anovan(iscSumAll,labelIndx, ...
%     'model', 'full', ...
%     'nest', nesting, ...
%     'random',[], ...
%     'varname', labels, 'display', 'on')
% saveas(figure(3),'stats','epsc')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% figure 6 - Component Comparisons
% clear h
% nModule = [1 2 3];
% indxViewing = [1 3 5; 2 4 6];
% nStim =size(indxViewing,2)
% plotIndx = repmat([1 1 2 2 3 3],1,nComp);
% barCol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]};
% iPlot = 1;
% i = 1
% 
% stimIndx = [1 2 3 ];
% visualIndx = [1 1 2];
% Allpower = []
% 
% f= figure(2);
% clf
% subplot 221
% for iComp = 1:3
%     for iStim = 1:nStim
%         attendPower = power{indxViewing(1,iStim)}(iComp,subjInt)';
%         disattendPower = power{indxViewing(2,iStim)}(iComp,:)';
%         
%         attendPowerMean = mean(attendPower,1);
%         disattendPowerMean = mean(disattendPower,1);
%         [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendPower,disattendPower);
%         % bar graph
%         hold on
%         
%         h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*(nStim*2+1),attendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
%             'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
%         iPlot = iPlot +1
%         h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*(nStim*2+1),disattendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
%             'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
%         
%         %graph property
%         yHeightMax = max([attendPowerMean disattendPowerMean]);
%         yHeight(i) = yHeightMax;
%         iPlot = iPlot +1
%         i = i+1;
%     end
% end
% 
% p =mystack(cell2mat(pValue'))
% pBonf = bonf_holm(p);
% pBonf(pBonf > .05) = nan
% pFdr =mafdr(p,'BHFDR',1);
% % pFdr(pFdr > .05) = nan;
% % p(p > .05) = nan
% % lineStyle = {'-k', '-k', '-k', '-k', '-k', '-k','-k', '-k', '-k'};
% % sigstarCustom({[1 2], [3 4], [5 6], [8 9], [10 11], [12 13], [15 16], [17 18],[19 20]},pFdr,[],lineStyle,yHeight+.001);
% 
% box Off
% set(gca, 'Xtick', [3.5 3.5+7 3.5+14], 'Xticklabel', {'C1', 'C2' ,'C3'});
% axis([0 21 .05 .125])
% l1=legend([h([1 3 5 7 9 11]) ], {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn' })
% ylabel('Alpha Power')
% rect = [0.525, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
% 
% figure(3)
% clf
% for i=1:6
%     attendPower(i,:) = power{i}(1,:);
% end
% 
% bar(attendPower');hold on
% 
% 
% for i=1:6
%     for ii=1:3
%     power_(ii,:,i)= power{i}(ii,:);
%     isc_(ii,:,i)=ISC{i}(ii,:);
%     end
%     rr(:,:,i)=corrcoef(power_(ii,:,i)',isc_(ii,:,i)');
% end

