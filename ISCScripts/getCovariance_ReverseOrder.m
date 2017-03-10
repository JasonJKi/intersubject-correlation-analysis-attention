% Get All Covariance.
clear all;
homeDir = '/home/jason/Repository/AttentionModulation/';
nModule = [1 2];

saveDir = 'Data/EEG/Processed/OrderReversed/'
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

groupIndx{1} = [1 2 3 4 5 6];
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
figure(1)
headprojection(A{1},3,'hor',0);
figure(2)
headprojection(A{2},3,'hor',0);


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
        ISC2{i} = isc;
        clear isc
            i= i+1;
    end
end

% Stim labels, and indexes for plotting
modality = [1 1 1 1 2 2];
attention = [1 2 1 2 1 2 ];
stimLabel = [1 1 2 2 3 3 ];
visualIndx = [1 1 1 1 0 0];
xAxis = [1 2 3 4 6 7 ];
narrative = [1 1 1 1 1 1];

% initialize variables
iscAll = []; iscSumAll = []; iscByCompAll =[]; iscRelativeContributionAll = [];
iscRelativeContributionMeanAll=[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
iscCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];

% set components of interest
compInt = 2 % components of interest
nComp = length(compInt) % measure first Ncomps 

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]};
figure(1); clf; sp1= subplot(5,5,[1 2 3 6 7 8 11 12 13]); 

stimIndx = [1 1 1 1 1 1 1];
subjgroup = [1 1 1 1 1 1 1];
stimOrder = [1 2 3 4 5 6];
nStim = length(stimOrder); deleteNegativeValues=0;
iPlot = 1;
subjInt=[1:3 6:12 15:20];
subjInt=[1:2 4 7:20];;

for iStim = 1:nStim
    isc_ = ISC2{stimOrder(iStim)}(compInt,subjInt);
    isc__ = ISC2{stimOrder(iStim)}(1:10,subjInt);
    
%     iscRelativeContribution = isc_./repmat(sum(abs(isc_),1),length(compInt),1);    
%     iscRelativeContributionAll = [iscRelativeContributionAll iscRelativeContribution];
%         
%     iscRelativeContributionMean = mean(iscRelativeContribution,2);
%     iscRelativeContributionMeanAll(:,iStim) = iscRelativeContributionMean;
%         
    iscSum = sum(isc_,1);
    iscByComp(:,iStim) = mean(isc__,2);
    iscSumMean(iStim) = mean(sum(isc_,1));
   
    if iStim == 1; iscSum(end) =[];end
    %labels
    nSubj=length(iscSum)
    subj = [1:nSubj]+((subjgroup(iStim))*size(iscSum,2)) ;
    attentionIndx = repmat(attention(iStim), nSubj,1);
    modalityIndx = repmat(modality(iStim), nSubj,1); 
    narrativeIndx = repmat(narrative(iStim),nSubj,1);
    
    if stimIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        iscSumAll = [iscSumAll; iscSum'];
        narrativeIndxAll = [narrativeIndxAll; narrativeIndx];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx]; 
        subjLabelAll = [subjLabelAll; subj'];
        % plottings
        hold on;
        if attention(iPlot) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),mean(iscSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        plot(xAxis(iPlot),iscSum,'.');
        n{iPlot} = num2str(length(iscSum));
        iPlot = iPlot +1;
    end
end

stimPair=[1 2; 3 4 ;5 6];
for i=1:3
    attentionDifference(i,:)=[sum(ISC2{stimPair(i,1)}(compInt,:),1) < sum(ISC2{stimPair(i,2)}(compInt,:),1)];
   [h(i) p(i)]=ttest([sum(ISC2{stimPair(i,1)}(compInt,:),1) - sum(ISC2{stimPair(i,2)}(compInt,:),1)]);

end

       
        

videos ={'BYD' 'GBU' 'PM'}
% plot properties
p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
l1=legend([p1 p2],'attend', 'count')
rect = [0.475, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF
delete([p1 p2]);
ylabel('ISC') ;
plotIndx = [1.5 3.5 6.5 8.5 10.5];
set(gca, 'Xtick', plotIndx,'xticklabel',videos)
t1=title('a');set(t1,'position',[-.05 1])


figure(2);clf
subplot(2,3,[1 2])
plotIndx = [1 1 2 2 3 3 1 1 2 2 3 3 1 1 2 2 3 3 ];
indxViewing=[1 3 5; 2 4 6];
barColIndx = [1 2 3]
barCol = [.3 .9 .9; 0 .6 .6;.5 .2 .8];
iPlot = 1;
i = 1

subjInt=[1:2 4 6:20];;

for iComp = 1:3
    for iStim = 1:3
        attendISC = ISC2{indxViewing(1,iStim)}(iComp,subjInt)';
        disattendISC = ISC2{indxViewing(2,iStim)}(iComp,subjInt)';
        attendISCMean = mean(attendISC,1);
        disattendISCMean = mean(disattendISC,1);
%         attendISCChanceMean = iscChanceMean(iComp,indxViewing(1,iStim)); 
%         disattendISCChanceMean = iscChanceMean(iComp,indxViewing(2,iStim));
         
        [sig(iComp,iStim) pValue{iComp,iStim}] = ttest(attendISC,disattendISC); 
        % bar graph
        hold on
        if (iComp == 3 & iStim == 3)
%             k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
%                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
            h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
            iPlot = iPlot +1
%             k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
%                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
            h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
        else
            h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
%             k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
%                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
            iPlot = iPlot +1
            h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean,0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
%             k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
%                 'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
        end
        
        %graph property
                
        yHeightMax = max([attendISCMean disattendISCMean]);
%         yHeightMax = max([attendISCMean disattendISCMean attendISCChanceMean disattendISCChanceMean]);
        yHeight(i) = yHeightMax;
        iPlot = iPlot +1
        i = i+1;
    end
end

p =mystack(cell2mat(pValue'))
pBonf = bonf_holm(p);
pBonf(pBonf > .05) = nan
pFdr =mafdr(p,'BHFDR',2);
pFdr(pFdr > .05) = nan;
p(p > .05) = nan
lineStyle = {'-k', '-k', '-k', '-k', '-k', '-k','-k', '-k', '-k'};
sigstarCustom({[1 2], [3 4], [5 6], [8 9], [10 11], [12 13], [15 16], [17 18],[19 20]},pBonf,[],lineStyle,yHeight+.001);
set(gca, 'Xtick', [3.5 10.5 17.5], 'Xticklabel', {'C1', 'C2' ,'C3'});
p1 =bar(20,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','-','LineWidth',1);
p2 =bar(21,1,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineStyle','--','LineWidth',1);
axis([0 21 0 .035])
l1=legend([h([2 4 6]) p2 p1], {'BYD' 'GBU' 'PM' 'count' 'attend'})
delete([p1 p2]);
ylabel('ISC')
box Off

saveas(figure(2),'Figure6','epsc')



%create alpha band filter
n = 7
fsRef =256
[b,a] = butter(n,[8 12]/fsRef*2,'bandpass');
freqz(b,a,fsRef,fsRef)
nComp = 3

% check the components
figure(1)
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
    for iSubj = 1:nSubj
        power_{i}(:,iSubj)= mean(abs(filter(b,a,eeg(:,:,iSubj)*W{compIndx(iModule)}(:,1:nComp)).^2),1);
        power_{i}(:,iSubj) =  power_{i}(:,iSubj)./mean(abs(eeg(:,:,iSubj)*W{compIndx(iModule)}(:,1:nComp)).^2,1)';
        disp(iSubj)
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


%% figure 6 - Component Comparisons
clear h
nModule = [1 2 3];
indxViewing = [1 3 5; 2 4 6];
nStim =size(indxViewing,2)
plotIndx = repmat([1 1 2 2 3 3],1,nComp);
barCol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]};
iPlot = 1;
i = 1

stimIndx = [1 2 3 ];
visualIndx = [1 1 2];
Allpower = []

f= figure(2);
clf
subplot 221
for iComp = 1:3
    for iStim = 1:nStim
        attendPower = power{indxViewing(1,iStim)}(iComp,subjInt)';
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
set(gca, 'Xtick', [3.5 3.5+7 3.5+14], 'Xticklabel', {'C1', 'C2' ,'C3'});
axis([0 21 .05 .125])
l1=legend([h([1 3 5 7 9 11]) ], {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn' })
ylabel('Alpha Power')
rect = [0.525, 0.85, .03, .06]; set(l1, 'Position', rect);legend BOXOFF

figure(3)
clf
for i=1:6
    attendPower(i,:) = power{i}(1,:);
end

bar(attendPower');hold on


for i=1:6
    for ii=1:3
    power_(ii,:,i)= power{i}(ii,:);
    isc_(ii,:,i)=ISC{i}(ii,:);
    end
    rr(:,:,i)=corrcoef(power_(ii,:,i)',isc_(ii,:,i)');
end