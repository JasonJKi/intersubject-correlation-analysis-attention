% compute the power at the Alpha Band.
clear all;
homeDir = '/home/jason/Repository/AttentionModulation/'; 
run([homeDir 'Codes/eegInfo.m'])
nModule = [1 2 3];
processtype = eeginfo.preprocess{2}
currentDir = 'Data/ISCValues/'
load([homeDir currentDir processtype(1:end-1) '_allComponents.mat'],'W','A')
compIndx = [1 1 2]; nComp = 10;

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
load([homeDir currentDir 'alpha_power_gabor2.mat'], 'power')

% Stim labels, and indexes for plotting
modality = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 ];
attention = [1 2 1 2 1 2 1 2 1 2 1 2 1 2];
stimLabel = [1 1 2 2 1 1 2 2 5 5 6 6 7 7];
visualIndx = [1 1 1 1 2 2 2 2 0 0 0 0];
xAxis = [1 2 3 4 6 7 8 9 11 12 13 14 15 16];
narrative = [1 1 1 1 1 1 1 1 1 1 2 2 2 2];

% initialize variables
powerAll = []; powerSumAll = []; powerByCompAll =[]; powerRelativeContributionAll = [];
powerRelativeContributionMeanAll=[]; 
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
powerCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];

% set components of interest
compInt = 1:3 % components of interest
nComp = length(compInt) % measure first Ncomps

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
stimIndx = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];
subjgroup = [1 1 1 1 2 2 2 2 3 3 3 3 3 3];
stimOrder = [1 3 2 4 5 7 6 8 9 12 10 13 11 14];
nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;
f =figure(1);
clf
subplot(2,9,[1:7])
for iStim = 1:nStim
    power_ = power{stimOrder(iStim)}(compInt,:);
    
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

videos = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
mode = {'audiovisual(free)' 'audiovisual(constrained)' 'audio only'} 
set(gca, 'Xtick',  [2.5 8 13.5 ], 'Xticklabel', mode);
% plot properties
p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
l1=legend([h([1 3 9 11 13]) p1 p2],[videos 'attend', 'count'])
rect = [0.76, 0.75, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
axis([min(xAxis)-1 max(xAxis)+1 .0 .6])
ylabel('alpha power')
box Off
saveas(f,'figure6a_attendVsDisattend','epsc')

%% compute significance figure 1
stimIndx = [1 2 5 6 9 10 11; 3 4 7 8 12 13 14];
nStim = size(stimIndx,2);
for iStim = 1:5
    attend = sum(power{stimIndx(1,iStim)}(compInt,:));
    disattend = sum(power{stimIndx(2,iStim)}(compInt,:));
    power_ = [attend disattend]';
    attentionIndx= [ones(length(attend),1) ;zeros(length(disattend),1)];
    [h1(iStim) p(iStim)] = ttest2(attend,disattend);
%     [p_{iStim} t]= anovan(power_,attentionIndx,'display','off');
    F(iStim) = t{2,6};
end

nesting = [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 1 0 0];
model = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 1 0; 0 0 0 1 1];
[P1,T1,STATS1,TERMS1] = anovan(powerSumAll,[subjLabelAll modalityIndxAll attentionIndxAll stimIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1 4], ...
    'varname', {'subject' 'modality' 'attention' 'stim'}, 'display', 'on')


%% figure 6 - Component Comparisons
clear h
nModule = [1 2 3];
indxViewing = [1 2 9 10 11; 3 4 12 13 14];
nStim =size(indxViewing,2)
plotIndx = repmat([1 1 2 2 3 3 4 4 5 5],1,nComp);
barCol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
iPlot = 1;
i = 1

stimIndx = [1 2 3 4 5];
visualIndx = [1 1 2 2 2];
Allpower = []

f= figure(2);
clf
subplot 221
for iComp = 1:3
    for iStim = 1:nStim
        attendPower = power{indxViewing(1,iStim)}(iComp,:)';
        disattendPower = power{indxViewing(2,iStim)}(iComp,:)';
        
        attendPowerMean = mean(attendPower,1);
        disattendPowerMean = mean(disattendPower,1);
        [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendPower,disattendPower);
        % bar graph
        hold on
        
        h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*(nStim*2+1),attendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
            'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
        iPlot = iPlot +1
        h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*(nStim*2+1),disattendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
            'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
        
        %graph property
        yHeightMax = max([attendPowerMean disattendPowerMean]);
        yHeight(i) = yHeightMax;
        iPlot = iPlot +1
        i = i+1;
    end
end

p =mystack(cell2mat(pValue'))
pBonf = bonf_holm(p);
pBonf(pBonf > .05) = nan
pFdr =mafdr(p,'BHFDR',1);
% pFdr(pFdr > .05) = nan;
% p(p > .05) = nan
% lineStyle = {'-k', '-k', '-k', '-k', '-k', '-k','-k', '-k', '-k'};
% sigstarCustom({[1 2], [3 4], [5 6], [8 9], [10 11], [12 13], [15 16], [17 18],[19 20]},pFdr,[],lineStyle,yHeight+.001);

box Off
set(gca, 'Xtick', [5.5 5.5+11 5.5+22], 'Xticklabel', {'C1', 'C2' ,'C3'});
axis([0 33 .05 .25])
l1=legend([h([1 3 5 7 9 11]) ], {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn' })
ylabel('Alpha Power')
rect = [0.525, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
saveas(f,'figure6b_attendVsDisattendByComponents','epsc')


%% compute significance for figure 6 b
compInt = 3; nComp = length(compInt);
narIndx = [1 1 1 2 2];
modalityIndx = [1 1 2 2 2];
stimIndxs = [1 2 3 4 5];
subjectGroupIndx = [1 1 2 2 2];

subjectIndxAllc = [];
attentionIndxAllc = [] ;
stimIndxAllc = [];
powerAllc= [];
componentAllc=[];
modalIndxAllc = [];

for iComp =1:3
    allpower = [];stimIndxAll = [];compIndxAll = [];
    attentionIndxAll = []; narrativeIndxAll =[];
    subjectIndxAll = []; modalIndxAll =[];
    for iStim = [1:5]
        attendpower = power{indxViewing(1,iStim)}(iComp,:)';
        disattendpower = power{indxViewing(2,iStim)}(iComp,:)';
        
        lAttend = size(attendpower,1);
        lDisattend = size(disattendpower,1);
        
        subjectIndx = [1:lAttend 1:lDisattend]'+(lAttend)*(subjectGroupIndx(iStim)-1);
        subjectIndxAll = [subjectIndxAll; subjectIndx];
        
        stimIndx = repmat(stimIndxs(iStim),(lAttend+lDisattend)*1,1);
        stimIndxAll = [stimIndxAll; stimIndx];
        
        modalIndx = repmat(modalityIndx(iStim),(lAttend+lDisattend)*1,1);
        modalIndxAll = [modalIndxAll; modalIndx];
        
        attentionIndx = [ones(lAttend*1,1); zeros(lDisattend*1,1)];
        attentionIndxAll = [attentionIndxAll; attentionIndx];
        
        power_ = [attendpower(:); disattendpower(:)];
        allpower = [allpower; power_];
    end

    anovaGroup =[subjectIndxAll modalIndxAll attentionIndxAll];
    nesting = [0 1 0; 0 0 0; 0 0 0];
    [P1{iComp},T1{iComp},STATS1{iComp},TERMS1{iComp}] = anovan(allpower,anovaGroup, ...
        'model', 'full', ...
        'nest', nesting, ...
        'random',[1], ...
        'varname', {'subject' 'modality' 'attention'}, ...
        'display', 'on')
    
    length(allpower)
    componentAllc = [componentAllc; repmat(iComp,length(allpower),1)]
    subjectIndxAllc = [subjectIndxAllc; subjectIndxAll];
    attentionIndxAllc = [attentionIndxAllc; attentionIndxAll];
    stimIndxAllc = [stimIndxAllc; stimIndxAll];
    powerAllc = [powerAllc; allpower];
    modalIndxAllc =[modalIndxAllc; modalIndxAll];
end

anovaGroup =[subjectIndxAllc modalIndxAllc attentionIndxAllc componentAllc stimIndxAllc];
nesting = [0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0];
[P1_,T1_,STATS1_,TERMS1_] = anovan(powerAllc,anovaGroup, ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1], ...
    'varname', {'subject' 'modality' 'attention' 'components' 'stimulus'}, ...
    'display', 'on')


for iComp =1:3
    allpower = [];stimIndxAll = [];compIndxAll = [];
    attentionIndxAll = []; narrativeIndxAll =[];
    subjectIndxAll = []; modalIndxAll =[];
    for iStim = [1:5]
        attendpower = power{indxViewing(1,iStim)}(iComp,:)';
        disattendpower = power{indxViewing(2,iStim)}(iComp,:)';
        
        lAttend = size(attendpower,1);
        lDisattend = size(disattendpower,1);
        
        subjectIndx = [1:lAttend 1:lDisattend]'+(lAttend)*(subjectGroupIndx(iStim)-1);
        subjectIndxAll = [subjectIndxAll; subjectIndx];
        
        stimIndx = repmat(stimIndxs(iStim),(lAttend+lDisattend)*1,1);
        stimIndxAll = [stimIndxAll; stimIndx];
        
        modalIndx = repmat(modalityIndx(iStim),(lAttend+lDisattend)*1,1);
        modalIndxAll = [modalIndxAll; modalIndx];
        
        attentionIndx = [ones(lAttend*1,1); zeros(lDisattend*1,1)];
        attentionIndxAll = [attentionIndxAll; attentionIndx];
        
        power_ = [attendpower(:); disattendpower(:)];
        allpower = [allpower; power_];
    end
    saveas(figure(2+iComp),['AlphaComponent' num2str(iComp) 'ANOVA'],'epsc')
    
    anovaGroup =[subjectIndxAll modalIndxAll attentionIndxAll];
    nesting = [0 1 0; 0 0 0; 0 0 0];
    [P1,T1,STATS1,TERMS1] = anovan(allpower,anovaGroup, ...
        'model', 'full', ...
        'nest', nesting, ...
        'random',[1], ...
        'varname', {'subject' 'modality' 'attention'}, ...
        'display', 'on')
end

%% figure 3 - Free Viewing vs Constrained Viewing
nModule = [1 2 3];

processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}

modality = [1 1 1 1 2 2 2 2 ];
attention = [1 2 1 2 1 2 1 2];
xAxis = [1 2 3 4 6 7 8 9 10];
narrative = [1 1 0 0 1 1 0 0];
subjgroup = [1 1 1 1 2 2 2 2];

power= power;
% initialize variables
powerAll = []; powerSumAll = []; powerByCompAll =[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
powerCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];
powerSumMean =[];
compInt = 1:3 % Components of Interest
nComp = length(compInt) % measure first Ncomps

% plot/ figure settings
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
stimLabel = [1 1 2 2 1 1 2 2];
stimOrder  = [1 3 2 4 5 7 6 8];
nStim = 8;
f = figure(3);clf
subplot 221
hold on
for iStim = 1:nStim
    power_ = power{stimOrder(iStim)}(compInt,:);
    powerSum = sum(power_,1);
    if iStim == 1; powerSum(end) =[];end
    powerSumMean(iStim) = mean(sum(power_,1));
    powerSumAll = [powerSumAll; powerSum']; %concat powers
    nSubj =  length(powerSum);
    nsubj_(iStim)=nSubj;
    % concat labels
    subj = [1:nSubj]+((subjgroup(iStim)-1)*size(powerSum,2)) ;
    subjLabelAll = [subjLabelAll; subj'];
    
    attentionIndx = repmat(attention(iStim), nSubj,1);
    attentionIndxAll = [attentionIndxAll; attentionIndx];
    
    stimIndx = repmat(stimLabel(iStim), nSubj,1);
    stimIndxAll = [stimIndxAll; stimIndx];
    
    modalityIndx = repmat(modality(iStim), nSubj,1);
    modalityIndxAll = [modalityIndxAll; modalityIndx];
    
    narrativeIndx = repmat(narrative(iStim),nSubj,1);
    narrativeIndxAll = [ narrativeIndxAll; narrativeIndx];
    % concat by components
    
    % plottings
    if attention(iStim) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
    h(iStim)= bar(xAxis(iStim),mean(powerSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
        'EdgeColor',barcol{stimLabel(iStim)}*barColorScale,'LineWidth',.75,'LineStyle',line,'BarWidth',.75);
    yheight1(iStim) = mean(powerSum);
    n_{iStim} = num2str(length(powerSum));
end

box Off
set(gca, 'Xtick', [2.5 7.5], 'Xticklabel', {'free viewing', 'constrained viewing'});
axis([0 10 .0 .4])
l1=legend([h([1 3]) ], {'BYD' 'GBU'})
ylabel('Alpha Power')
rect = [0.50, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
saveas(f,'figure6c_attendVsDisattendVisualConstraints','epsc')

figure(8)
nesting = [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 1 0 0];
model = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 1 0; 0 0 0 1 1];
[P1,T1,STATS1,TERMS1] = anovan(powerSumAll,[subjLabelAll modalityIndxAll stimIndxAll attentionIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1], ...
    'varname', {'subjects' 'visual condition' 'attention' 'stimulus'}, 'display', 'on');

saveas(figure(4),'AlphafreeVsfixationANOVA','epsc')



%% figure 2 - Narrative vs Non-Narrative Component Comparisons
nModule = [1 2 3];
indxViewing = [1 2 5 6 9 10 11; 3 4 7 8 12 13 14];
figure(2);clf
subplot(2,3,[1 2])
nStim =size(indxViewing,2)
plotIndx = repmat([1 1 2 2 1 1 2 2 3 3 4 4 5 5],1,nComp);
barCol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
iPlot = 1;
i = 1


stimIndx = [1 2 3 4 5];
visualIndx = [1 1 2 2 2];
Allpower = []
for iComp = 1:3
    for iStim = 1:nStim
        attendPower = power{indxViewing(1,iStim)}(iComp,:)';
        disattendPower = power{indxViewing(2,iStim)}(iComp,:)';
        
        attendPowerMean = mean(attendPower,1);
        disattendPowerMean = mean(disattendPower,1);
        [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendPower,disattendPower);
        % bar graph
        hold on
        
        h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*(nStim*2+1),attendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
            'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
        iPlot = iPlot +1
        h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*(nStim*2+1),disattendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
            'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
        
        %graph property
        yHeightMax = max([attendPowerMean disattendPowerMean]);
        yHeight(i) = yHeightMax;
        iPlot = iPlot +1
        i = i+1;
    end
end
stimOrder = [1 3 2 4 5 7 6 8 9 12 10 13 11 14];
nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;
figure(6);
for iStim = 1:nStim
    power_ = power{stimOrder(iStim)}(compInt,:);
    
    powerSum = sum(power_,1);
    powerSumMean(iStim) = mean(sum(power_,1));
    
    if iStim == 1; powerSum(end) =[];end
    %labels
    nSubj=length(powerSum)
    subj = [1:nSubj]+((subjgroup(iStim))*size(powerSum,2)) ;
    attentionIndx = repmat(attention(iStim), nSubj,1);
    modalityIndx = repmat(modality(iStim), nSubj,1);
    narrativeIndx = repmat(narrative(iStim),nSubj,1);
    
    if stimIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        powerSumAll = [powerSumAll; powerSum'];
        narrativeIndxAll = [narrativeIndxAll; narrativeIndx];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx];
        subjLabelAll = [subjLabelAll; subj'];
        % plottings
        hold on;
        if attention(iPlot) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),mean(powerSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        %         plot(xAxis(iPlot),powerSum,'.');
        iPlot = iPlot +1
    end
end

saveas(figure(6),'alphaModality','epsc')

%% figure 2 - Component Comparisons

load([homeDir currentDir 'alpha_power_butterworth2_unnormalized.mat'], 'power')
clear h
nModule = [1 2 3];
indxViewing = [1 2 5 6 9 10 11; 3 4 7 8 12 13 14];
nStim =size(indxViewing,2)
plotIndx = repmat([1 1 2 2 3 3 4 4 5 5 6 6 7 7],1,nComp);
barCol = {[.3 .9 .9];[ 0 .6 .6]; [.4 .8 .8];[ .1 .7 .7]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
iPlot = 1;
i = 1

stimIndx = [1 2 3 4 5];
visualIndx = [1 1 2 2 2];
Allpower = []

f= figure(2);
clf
subplot 221
for iComp = 1:3
    for iStim = 1:nStim
        attendPower = power{indxViewing(1,iStim)}(iComp,:)';
        disattendPower = power{indxViewing(2,iStim)}(iComp,:)';
        
        attendPowerMean = mean(attendPower,1);
        disattendPowerMean = mean(disattendPower,1);
        [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendPower,disattendPower);
        % bar graph
        hold on
        
        h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*(nStim*2+1),attendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
            'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
        iPlot = iPlot +1
        h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*(nStim*2+1),disattendPowerMean,0.5,'FaceColor', barCol{plotIndx(iPlot),:}, ...
            'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
        
        %graph property
        yHeightMax = max([attendPowerMean disattendPowerMean]);
        yHeight(i) = yHeightMax;
        iPlot = iPlot +1
        i = i+1;
    end
end

p =mystack(cell2mat(pValue'))
pBonf = bonf_holm(p);
pBonf(pBonf > .05) = nan
pFdr =mafdr(p,'BHFDR',1);
% pFdr(pFdr > .05) = nan;
% p(p > .05) = nan
% lineStyle = {'-k', '-k', '-k', '-k', '-k', '-k','-k', '-k', '-k'};
% sigstarCustom({[1 2], [3 4], [5 6], [8 9], [10 11], [12 13], [15 16], [17 18],[19 20]},pFdr,[],lineStyle,yHeight+.001);

box Off
set(gca, 'Xtick', [7.5 7.5+15 7.5+30], 'Xticklabel', {'C1', 'C2' ,'C3'});
axis([0 45 .00 .10])
l1=legend([h([1 3 5 7 9 11 13 15]) ], {'BYD free' 'GBU free' 'BYD const.' 'GBU const.' 'PM' 'PMsc' 'Jpn' })
ylabel('Alpha Power')
rect = [0.525, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
saveas(f,'figure6b_attendVsDisattendByComponentsAll_unnormalized','epsc')

%% figure 3 - power-A and Attentional Prediction)

%Generate figure 3
figure(3)
% clf
% compInt = 1:3;
% bng_1 = sum(power{1}(compInt,:),1)';
% bng_2 = sum(power{3}(compInt,:),1)';
% gbu_1 = sum(power{2}(compInt,:),1)';
% gbu_2 = sum(power{4}(compInt,:),1)';
% 
% for j = 1:size(bng_1,2);
%     text(.8,bng_1(j),num2str(j));
% end
% 
% (bng_2'-bng_1' > 0)
% (bng_2'-bng_1' > 0)
% 
% [Azbng,tpbng,fpbng,fcbng] = rocarea([bng_2' bng_1'], [ones(1,length(bng_1)) zeros(1,length(bng_2))]);
% [Azgbu,tpgbu,fpgbu,fcgbu] = rocarea([gbu_2' gbu_1' ], [ones(1,length(gbu_1)) zeros(1,length(gbu_2))]);
% azgbu = num2str(Azgbu);
% azbng = num2str(Azbng);
% 
% subplot(2,2,1)
% hold all
% h2= plot([1 2], [bng_1 bng_2], 'Color' ,[.3 .9 .9],'Marker' , 'x');
% h3= plot([3 4], [gbu_1 gbu_2], 'Color' , [0 .6 .6], 'Marker' , 'x');
% hold off
% k=legend([h2(1), h3(1)],[' BYD (Az = ' azbng(1:end) ')'],[' GBU (Az = ' azgbu  ')'])
% legend BOXOFF
% rect = [0.35, 0.85, .025, .025];
% set(k, 'Position', rect)
% 
% ylabel('power');
% XMIN = .5;XMAX = 4.5;YMIN = 0;YMAX = max(gbu_2)+.01;
% axis([XMIN XMAX YMIN YMAX])
% set(gca, 'Xtick', [1 2 3 4] ,'xticklabel', {'attend' 'count' 'attend' 'count'});
% box OFF

attention =[1 2 9 10 11; 3 4 12 13 14];
stimOrder= [1 3 6 8 10; 2 4 7 9 11];
azIndx = [1 2 4 5 6];

subplot(2,2,1)
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
%Az Significance.
Nshuffles = 10000;
for iStim =1:5
    powerAttend = sum(power{attention(1,iStim)}(compInt,:));
    powerCount = sum(power{attention(2,iStim)}(compInt,:));
    label = [ones(length(powerAttend),1); zeros(length(powerCount),1)];
    az(iStim) = rocarea([powerAttend powerCount],label);
    lengthlabel = length(label);
    for iShuffle = 1:Nshuffles
        az_(iShuffle) = rocarea([powerAttend powerCount], label(randperm(lengthlabel)));
    end
    p = length(find(az_ > az(iStim)))/Nshuffles;
    i = 10
    rnd = 1
    while rnd
        if p == 0
            ps = num2str(.0001,'%.4f');
            rnd =0
        elseif p < 1/10^(i)
            i
            p_=floor(p*10^(i+2))/10^(i+2);
            ps= num2str(p_);
            rnd = 0;
        end
        i= i-1
    end
    P{iStim} = num2str(ps,'%.4f');
    cdfAz = cdf(makedist('Normal',mean(az_),std(az_)),az_);
    sigPvalue(iStim) = min(az_(find(cdfAz> .95)));
    clear az_
end
sigPvalue = mean(sigPvalue);

% shade area of insignificance
x=[0 (nStim+2)]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
shadedplot([0 7], y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on

for iStim =1:5
    hold on
    powerAttend = sum(power{attention(2,iStim)}(compInt,:),1);
    powerCount = sum(power{attention(1,iStim)}(compInt,:),1);
    label = [ones(length(powerAttend),1); zeros(length(powerCount),1)];
    az(iStim) = rocarea([powerAttend powerCount], label);
    bar(azIndx(iStim), az(iStim),0.5,'FaceColor', barcol{iStim}, ...
         'EdgeColor',barcol{iStim}*.02,'LineWidth',1,'BarWidth',.75)
end

plot([0 7],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1,'Color', [.25 .25 .25])
ylabel('Az');stims = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
set(gca, 'Xtick', azIndx ,'xticklabel', stims, 'Ytick', [.5 .75 1], 'Yticklabel', {'0.5' '.75' '1'});
axis([0 7 .5 1])
saveas(figure(3), 'figure6b', 'epsc')

%% normalize subjects
figure(4)
clear subjAll
attention =[1 3 5 7 9; 2 4 6 8 10];
stimIndx = [1 3 2 4 9 12 10 13 11 14];
nStim = length(stimIndx)
compInt = 1:3;

for i = 1:4
    subjAll(:,i) = sum(power{stimIndx(i)}(compInt,:),1)';
end
subjMean = mean(subjAll,2);
subjAll = subjAll./repmat(subjMean,1,4);

for i = 5:10
    subjAll(:,i) = sum(power{stimIndx(i)}(compInt,:),1)';
end
subjMean = mean(subjAll(:,5:10),2);
subjAll(:,5:10) = subjAll(:,5:10)./repmat(subjMean,1,6);

[Azbng,tpbng,fpbng,fcbng] = rocarea([subjAll(:,2);subjAll(:,1)], [ones(1,length(subjAll(:,1))) zeros(1,length(subjAll(:,2)))]);
[Azgbu,tpgbu,fpgbu,fcgbu] = rocarea([subjAll(:,4);subjAll(:,3)], [ones(1,length(subjAll(:,3))) zeros(1,length(subjAll(:,4)))]);
azgbu = num2str(Azgbu);
azbng = num2str(Azbng);

subplot(2,2,1)
% hold all
% h2= plot([1 2], [subjAll(:,1) subjAll(:,2)], 'Color' ,[.3 .9 .9],'Marker' , 'x');
% h3= plot([3 4], [subjAll(:,3) subjAll(:,4)], 'Color' , [0 .6 .6], 'Marker' , 'x');
% hold off
% k=legend([h2(1), h3(1)],[' BYD (Az = ' azbng(1:end) ')'],[' GBU (Az = ' azgbu  ')'])
% legend BOXOFF
% rect = [0.35, 0.85, .025, .025];
% set(k, 'Position', rect)
% ylabel('power');
% XMIN = .5;XMAX = 4.5;YMIN = 0;YMAX = max(subjAll(:,4))+.01;
% axis([XMIN XMAX YMIN YMAX])
% set(gca, 'Xtick', [1 2 3 4] ,'xticklabel', {'attend' 'count' 'attend' 'count'});box OFF;
% 
% subplot(2,2,2)
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
%Az Significance.
Nshuffles = 10000;
for iStim =1:5
    powerAttend = subjAll(:,attention(1,iStim));
    powerCount = subjAll(:,attention(2,iStim));
    label = [ones(length(powerAttend),1); zeros(length(powerCount),1)];
    az(iStim) = rocarea([powerCount;powerAttend],label);
    lengthlabel = length(label);
    for iShuffle = 1:Nshuffles
        az_(iShuffle) = rocarea([powerAttend powerCount], label(randperm(lengthlabel)));
    end
    p = length(find(az_ > az(iStim)))/Nshuffles;
    i = 10
    rnd = 1
    while rnd
        if p == 0
            ps = num2str(.0001,'%.4f');
            rnd =0
        elseif p < 1/10^(i)
            i
            p_=floor(p*10^(i+2))/10^(i+2);
            ps= num2str(p_);
            rnd = 0;
        end
        i= i-1
    end
    P{iStim} = num2str(ps,'%.4f');
    cdfAz = cdf(makedist('Normal',mean(az_),std(az_)),az_);
    sigPvalue(iStim) = min(az_(find(cdfAz> .95)));
    clear az_
end

sigPvalue = mean(sigPvalue);
% shade area of insignificance
x=[0 (nStim+2)]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
shadedplot([0 7], y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on

for iStim =1:5
    hold on
   powerAttend = subjAll(:,attention(1,iStim));
    powerCount = subjAll(:,attention(2,iStim));
     label = [ones(length(powerAttend),1); zeros(length(powerCount),1)];
    az(iStim) = rocarea([powerCount; powerAttend], label);
    bar(azIndx(iStim), az(iStim),0.5,'FaceColor', barcol{iStim}, ...
         'EdgeColor',barcol{iStim}*.02,'LineWidth',1,'BarWidth',.75)
end

plot([0 7],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1,'Color', [.25 .25 .25])
ylabel('Az');stims = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
set(gca, 'Xtick', azIndx ,'xticklabel', stims, 'Ytick', [.5 .75 1], 'Yticklabel', {'0.5' '.75' '1'});
axis([0 7 .5 1])
saveas(figure(4), 'figure6b_normalized', 'epsc')

%%

clear subjAll
attention =[1 3 5 7 9 11 13 15; 2 4 6 8 10 12 14 16];
stimIndx = [1 3 2 4 5 7 6 8 9 12 10 13 11 14];
nStim = length(stimIndx)
compInt = 1:3;

for i = 1:4
    subjAll(:,i) = sum(power{i}(compInt,:),1)';
end
subjMean = mean(subjAll,2);
subjAll = subjAll./repmat(subjMean,1,4);

for i = 1:4
    powerNormalized{i} = subjAll(:,i);
end
clear subjAll

ii = 1    
for i = 5:8
    subjAll(:,ii) = sum(power{i}(compInt,:),1)';
    ii = ii+1;
end
subjMean = mean(subjAll,2);
subjAll = subjAll./repmat(subjMean,1,4);

ii = 1
for i = 5:8
    powerNormalized{i} = subjAll(:,ii);
    ii = ii+1;
end
clear subjAll

ii = 1
for i = 9:14
    subjAll(:,ii) = sum(power{i}(compInt,:),1)';
    ii = ii+1;
end
subjMean = mean(subjAll,2);
subjAll = subjAll./repmat(subjMean,1,6);


ii = 1
for i = 9:14
    powerNormalized{i} = subjAll(:,ii);
    ii = ii+1;
end
clear subjAll

% Stim labels, and indexes for plotting
modality = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 ];
attention = [1 2 1 2 1 2 1 2 1 2 1 2 1 2];
stimLabel = [1 1 2 2 3 3 4 4 5 5 6 6 7 7];
visualIndx = [1 1 1 1 2 2 2 2 0 0 0 0];
xAxis = [1 2 3 4 6 7 8 9 11 12 13 14 15 16];
narrative = [1 1 1 1 1 1 1 1 1 1 2 2 2 2];

% initialize variables
powerAll = []; powerSumAll = []; powerByCompAll =[]; powerRelativeContributionAll = [];
powerRelativeContributionMeanAll=[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
powerCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];

% set components of interest
compInt = 1:3 % components of interest
nComp = length(compInt) % measure first Ncomps

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
stimIndx = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];
subjgroup = [1 1 1 1 2 2 2 2 3 3 3 3 3 3];
stimOrder = [1 3 2 4 5 7 6 8 9 12 10 13 11 14];
nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;
f =figure(1);
clf
subplot(2,9,[1:7])
for iStim = 1:nStim

    powerSum = mean(powerNormalized{stimOrder(iStim)});
    powerSumMean(iStim) = mean(powerSum,1);
    
%     if iStim == 1; powerSum(end) =[];end
    %labels
    nSubj=length(powerSum)
    subj = [1:nSubj]+((subjgroup(iStim))*size(powerSum,2)) ;
    attentionIndx = repmat(attention(iStim), nSubj,1);
    modalityIndx = repmat(modality(iStim), nSubj,1);
    narrativeIndx = repmat(narrative(iStim),nSubj,1);
    
    if stimIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        powerSumAll = [powerSumAll; powerSum'];
        narrativeIndxAll = [narrativeIndxAll; narrativeIndx];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx];
        subjLabelAll = [subjLabelAll; subj'];
        % plottings
        hold on;
        if attention(iPlot) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),powerSum,0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        %                 plot(xAxis(iPlot),powerSum,'.');
        iPlot = iPlot +1
    end
end

videos = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
mode = {'audiovisual(free)' 'audiovisual(constrained)' 'audio only'} 
set(gca, 'Xtick',  [2.5 8 13.5 ], 'Xticklabel', mode);
% plot properties
p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
l1=legend([h([1 3 9 11 13]) p1 p2],[videos 'attend', 'count'])
rect = [0.76, 0.75, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
axis([min(xAxis)-1 max(xAxis)+1 .3 .75])
ylabel('alpha power')
box Off
