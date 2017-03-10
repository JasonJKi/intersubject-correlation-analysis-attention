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
clear all; close all
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 1 - ATTENTION COMPARISON (BYD,GBU,PM,PMsc,Jpn)
% Stim labels, and indexes for plotting
modality = [1 1 1 1 1 1 1 1 2 2 2 2 2 2 ];
attention = [1 2 1 2 1 2 1 2 1 2 1 2 1 2];
stimLabel = [1 1 2 2 1 1  2 2 3 3 4 4 5 5];
visualIndx = [1 1 1 1 2 2 2 2 0 0 0 0];
xAxis = [1 2 3 4 6 7 8 9 10 11];
narrative = [1 1 1 1 1 1 1 1 1 1 2 2 2 2];

% initialize variables
iscAll = []; iscSumAll = []; iscByCompAll =[]; iscRelativeContributionAll = [];
iscRelativeContributionMeanAll=[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
iscCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];stimulusIndxAll=[];

% set components of interest
compInt = 1:3 % components of interest
nComp = length(compInt) % measure first Ncomps 

% figure and plotting properties
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
figure(1); clf; sp1= subplot(5,5,[1 2 3 6 7 8 11 12 13]); 

stimIndx = [1 1 1 1 0 0 0 0 1 1 1 1 1 1];
subjgroup = [1 1 1 1 0 0 0 0 2 2 2 2 2 2];
stimOrder = [1 3 2 4 5 7 6 8 9 12 10 13 11 14];
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
    narrativeIndx = repmat(narrative(iStim),nSubj,1);
    stimulusIndx = repmat(stimLabel(iStim),nSubj,1);

    if stimIndx(iStim) %plot stimulus of interest
        %concat labels of interest
        iscSumAll = [iscSumAll; iscSum'];
        narrativeIndxAll = [narrativeIndxAll; narrativeIndx];
        modalityIndxAll = [modalityIndxAll; modalityIndx];
        attentionIndxAll = [attentionIndxAll; attentionIndx]; 
        stimulusIndxAll= [stimulusIndxAll; stimulusIndx];
        subjLabelAll = [subjLabelAll; subj'];
        % plottings
        hold on;
        if attention(iPlot) == 1; barColorScale =1/2;line= '-';else barColorScale = 1/4;line= '--'; end
        h(iPlot)= bar(xAxis(iPlot),mean(iscSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
            'EdgeColor',[0 0 0],'LineWidth',1.3,'LineStyle',line,'BarWidth',.75);
        plot(xAxis(iPlot),iscSum,'.');
        n{iPlot} = num2str(length(iscSum));
        iPlot = iPlot +1
    end
end

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

% compute significance 
stimIndx = [1 2 9 10 11; 3 4 12 13 14];
for iStim = 1:5
    attend = sum(ISC{stimIndx(1,iStim)}(compInt,:));
    disattend = sum(ISC{stimIndx(2,iStim)}(compInt,:));
    isc = [attend disattend]';
    attentionIndx= [ones(length(attend),1) ;zeros(length(disattend),1)];
    [h1(iStim) p(iStim) CI{iStim} STATS{iStim}] = ttest(attend,disattend);
    [p_{iStim} t]= anovan(isc,attentionIndx,'display','off');
    F(iStim) = t{2,6};
end

linestyle = {'-k' '-k' '-k' '-k' '-k'};
yheight = [0 0 .065 0 0];
pBonf = bonf_holm(p);
pBonf(pBonf > .05) = nan;
pFdr =mafdr(p,'BHFDR',1);
pFdr(pFdr > .05) = nan;
p(p > .05) = nan
sigstarCustom({[1,2],[3,4], [6,7], [8, 9], [10,11]},pBonf,[],linestyle,yheight)
saveas(figure(1),[homeDir 'Data/figure1'],'epsc')

%% Compute statistical significance of the components
load('iscChanceAll')
iscChanceAll = permute(iscChanceAll,[2 1 3]); whos iscChanceAll
for i = 1:14
    for ii = 1:10
    mu(i,ii) = mean(iscChanceAll(:,ii,stimOrder(i)));
    sigma = std(iscChanceAll(:,ii,stimOrder(i)));
    y(ii,i) = cdf('Normal',iscByComp(ii,i),mu(i,ii),sigma);
    end
end

for i = 1:14
    for ii = 1:10
        iscChanceMean(ii,i) = mean(iscChanceAll(:,ii,i));
        pVal(ii,i) =  sum(iscChanceAll(:,ii,i) > iscByComp(ii,i))/size(iscChanceAll,1)
    end
end
% save('iscChanceMean')

nesting = [0 1 0 0; 0 0 0 0; 0 0 0 0;0 1 0 0];
% model = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 1 0; 0 0 0 1 1];
[P1,T1,STATS1,TERMS1] = anovan(iscSumAll,[subjLabelAll modalityIndxAll attentionIndxAll stimulusIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1], ...
    'varname', {'subject' 'modality' 'attention' 'stimulus'}, 'display', 'on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 2 - Component Comparisons
homeDir = '/home/jason/Repository/AttentionModulation/';
run([homeDir 'Codes/eegInfo.m'])
nModule = [1 2 3];

processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
iscFile =  [processType(1:end-1) '_allISC.mat'];
saveDir = 'Data/ISCValues/'
saveDir = [homeDir saveDir];
load([saveDir iscFile]);
indxViewing = [9 10 11;12 13 14];

xx=load('iscChanceMean')
iscChanceMean=xx.iscChanceMean;

figure(2);clf
subplot(2,3,[1 2])
plotIndx = [1 1 2 2 3 3 1 1 2 2 3 3 1 1 2 2 3 3 ];
barColIndx = [1 2 3]
barCol = [.5 .2 .8; .6 .4 .5; .7 .6 .5];
iPlot = 1;
i = 1
for iComp = 1:3
    for iStim = 1:3
        attendISC = ISC{indxViewing(1,iStim)}(iComp,:)';
        disattendISC = ISC{indxViewing(2,iStim)}(iComp,:)';
        attendISCMean(iComp,iStim) = mean(attendISC,1);
        disattendISCMean(iComp,iStim) = mean(disattendISC,1);
        attendISCChanceMean = iscChanceMean(iComp,indxViewing(1,iStim)); 
        disattendISCChanceMean = iscChanceMean(iComp,indxViewing(2,iStim));
        
        [sig(iComp,iStim) pValue{iComp,iStim} CI{iComp,iStim} STATS{iComp,iStim}] = ttest(attendISC,disattendISC); 
        % bar graph
        hold on
        if (iComp == 3 & iStim == 3)
            k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
                'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
            h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean(iComp,iStim),0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
            iPlot = iPlot +1
            k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
                'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
            h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean(iComp,iStim),0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
        else
            h(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCMean(iComp,iStim),0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','-','BarWidth',.75);
            k(iPlot) = bar(1+(iStim-1)*2+(iComp-1)*7,attendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
                'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
            iPlot = iPlot +1
            h(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCMean(iComp,iStim),0.5,'FaceColor', barCol(plotIndx(iPlot),:), ...
                'EdgeColor',[0 0 0],'LineWidth',1,'LineStyle','--','BarWidth',.75);
            k(iPlot)= bar(2+(iStim-1)*2+(iComp-1)*7,disattendISCChanceMean,0.5,'FaceColor', [.2 .2 .2], ...
                'EdgeColor',[0 0 0],'LineWidth',.1,'LineStyle','-','BarWidth',.75);
        end
        
        %graph property
        yHeightMax = max([attendISCMean(iComp,iStim) disattendISCMean(iComp,iStim) attendISCChanceMean disattendISCChanceMean]);
        yHeight(i) = yHeightMax;
        iPlot = iPlot +1
        i = i+1;
    end
end

p =mystack(cell2mat(pValue'))
pfdr = bonf_holm(p);
pfdr =mafdr(p,'BHFDR',2);
pfdr(pfdr > .05) = nan;
p(p > .05) = nan
lineStyle = {'-k', '-k', '-k', '-k', '-k', '-k','-k', '-k', '-k'};
sigstarCustom({[1 2], [3 4], [5 6], [8 9], [10 11], [12 13], [15 16], [17 18],[19 20]},pfdr,[],lineStyle,yHeight+.001);
set(gca, 'Xtick', [3.5 10.5 17.5], 'Xticklabel', {'C1', 'C2' ,'C3'});
axis([0 21 0 .025])
l1=legend([h([1 3 5]) ], {'PM' 'PMsc' 'Jpn' })
ylabel('ISC')
box Off

saveas(figure(2),[homeDir 'Data/figure2'],'epsc');


%% significance
compInt = 1;nComp = length(compInt);
narIndx = [1 2 2];

allISC = [];stimIndxAll = [];compIndxAll = [];
attentionIndxAll = []; narrativeIndxAll =[];
subjectIndxAll = [];
for iStim = 1:3
    attendISC = ISC{indxViewing(1,iStim)}(compInt,:)';
    disattendISC = ISC{indxViewing(2,iStim)}(compInt,:)';
    
    lAttend = size(attendISC,1);
    lDisattend = size(disattendISC,1);
    
    subjectIndx = [repmat(1:lAttend,1,3) repmat(1:lDisattend,1,3)]';
    subjectIndxAll = [subjectIndxAll; subjectIndx];
    
    stimIndx = repmat(iStim,(lAttend+lDisattend)*nComp,1);
    stimIndxAll = [stimIndxAll; stimIndx];
    
    compIndxA = repmat([1 2 3],lAttend,1); 
    compIndxB = repmat([1 2 3],lDisattend,1);
    compIndx = [compIndxA(:); compIndxB(:)];
    compIndxAll = [compIndxAll; compIndx];
    
    narrativeIndx = repmat(narIndx(iStim),(lAttend+lDisattend)*nComp,1);
    narrativeIndxAll = [narrativeIndxAll; narrativeIndx];
    
    attentionIndx = [ones(lAttend*nComp,1); zeros(lDisattend*nComp,1)];
    attentionIndxAll = [attentionIndxAll; attentionIndx];
    
    isc = [attendISC(:); disattendISC(:)];
    allISC = [allISC; isc];
end

anovaGroup =[subjectIndxAll stimIndxAll attentionIndxAll];
nesting = [0 0 0; 0 0 0; 0 0 0];
model = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 1 0; 0 0 0 1 1];
[P1,T1,STATS1,TERMS1] = anovan(allISC,anovaGroup, ...
    'model', 'full', ...
    'nest', [], ...
    'random',[1], ...
    'varname', {'subject' 'stimulus' 'attention'}, ...
    'display', 'on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 3 - ISC-A and Attentional Prediction)
processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
iscFile =  [processType(1:end-1) '_allISCPredictor.mat'];
saveDir = 'Data/ISCValues/'
saveDir = [homeDir saveDir];
load([saveDir iscFile]);

%Generate figure 3
figure(3)
clf
compInt = 1:3;
bng_1 = sum(ISC{1}(compInt,:),1)';
bng_2 = sum(ISC{3}(compInt,:),1)';
gbu_1 = sum(ISC{2}(compInt,:),1)';
gbu_2 = sum(ISC{4}(compInt,:),1)';

for j = 1:size(bng_1,2);
    text(.8,bng_1(j),num2str(j));
end

[Azbng,tpbng,fpbng,fcbng] = rocarea([bng_1' bng_2'], [ones(1,length(bng_1)) zeros(1,length(bng_2))]);
[Azgbu,tpgbu,fpgbu,fcgbu] = rocarea([gbu_1' gbu_2'], [ones(1,length(gbu_1)) zeros(1,length(gbu_2))]);
azgbu = num2str(Azgbu);
azbng = num2str(Azbng);

subplot(2,2,1)
hold all
h2= plot([1 2], [bng_1 bng_2], 'Color' ,[.3 .9 .9],'Marker' , 'x');
h3= plot([3 4], [gbu_1 gbu_2], 'Color' , [0 .6 .6], 'Marker' , 'x');
hold off
k=legend([h2(1), h3(1)],[' BYD (Az = ' azbng(1:end) ')'],[' GBU (Az = ' azgbu  ')'])
legend BOXOFF
rect = [0.35, 0.85, .025, .025];
set(k, 'Position', rect)

ylabel('ISC-A');
XMIN = .5;XMAX = 4.5;YMIN = 0;YMAX = max(bng_1)+.01;
axis([XMIN XMAX YMIN YMAX])
set(gca, 'Xtick', [1 2 3 4] ,'xticklabel', {'attend' 'count' 'attend' 'count'});
box OFF

attention =[1 2 9 10 11; 3 4 12 13 14];
stimOrder= [1 3 6 8 10; 2 4 7 9 11];
azIndx = [1 2 4 5 6];

%Az Significance.
Nshuffles = 10000;
for iStim =1:5
    iscAttend = sum(ISC{attention(1,iStim)}(compInt,:));
    iscCount = sum(ISC{attention(2,iStim)}(compInt,:));
    label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
    az(iStim) = rocarea([iscAttend iscCount],label);
    lengthlabel = length(label);
    for iShuffle = 1:Nshuffles
        az_(iShuffle) = rocarea([iscAttend iscCount], label(randperm(lengthlabel)));
    end
    p = length(find(az_ > az(iStim)))/Nshuffles;
    pp(iStim)=p
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
end

cdfAz = cdf(makedist('Normal',mean(az_),std(az_)),az_)
sigPvalue = min(az_(find(cdfAz> .99))) 

% shade area of insignificance
subplot(2,2,2)
x=[0 (nStim+2)]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
shadedplot([0 7], y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};

for iStim =1:5
    hold on
    iscAttend = sum(ISC{attention(1,iStim)}(compInt,:),1);
    iscCount = sum(ISC{attention(2,iStim)}(compInt,:),1);
    label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
    az(iStim) = rocarea([iscAttend iscCount], label);
    bar(azIndx(iStim), az(iStim),0.5,'FaceColor', barcol{iStim}, ...
         'EdgeColor',barcol{iStim}*.02,'LineWidth',1,'BarWidth',.75)
end

plot([0 7],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1,'Color', [.25 .25 .25])
ylabel('Az');stims = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
set(gca, 'Xtick', azIndx ,'xticklabel', stims, 'Ytick', [.5 .75 1], 'Yticklabel', {'0.5' '.75' '1'});
axis([0 7 .5 1])
saveas(figure(3), 'figure3', 'epsc')

%% figure4 - electrode position Analysis
figure(4)
processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
iscFile =  [processType(1:end-1) '_allISCPredictor.mat'];
iscFile = ['allISCAnteroposteriorAxis'];
saveDir = 'Data/ISCValues/'
saveDir = [homeDir saveDir];


load([saveDir iscFile]);
nElectrodes= 9;nComp = 1:3;sumISCFree = [];sumISCCount = [];AzAll_=[];

linewidth = [2 2 1 1 ]
attentionIndx = [1 2 9;3 4 12];
nStim = size(attentionIndx,2);
for iElectrode = 1:nElectrodes
    for iStim = 1:nStim
        ISCFree = sum(ISC{attentionIndx(1,iStim)}{iElectrode}(nComp,:),1);
        ISCCount =sum(ISC{attentionIndx(2,iStim)}{iElectrode}(nComp,:),1);
        Az(iStim) = rocarea([ISCFree ISCCount],[ones(1,length(ISCFree)) zeros(1,length(ISCCount))]);
    end
    AzAll(iElectrode,:) = Az;
    clear Az
end

% compute significance
iStim = 1;compInt =1:3;iElectrode = 1; Nshuffles = 10000
iscAttend = sum(ISC{attentionIndx(1,iStim)}{iElectrode}(compInt,:));
iscCount = sum(ISC{attentionIndx(2,iStim)}{iElectrode}(compInt,:));
label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
lengthlabel = length(label);
for iShuffle = 1:Nshuffles
    az_(iShuffle) = rocarea([iscAttend iscCount], label(randperm(lengthlabel)));
end

cdfAz = cdf(makedist('Normal',mean(az_),std(az_)),az_);
sigPvalue = min(az_(find(cdfAz> .99)));

%color coding and plots
barcol = {[.3 .9 .9]; [0 .6 .6];  [.5 .2 .8]};

figure(4)
clf
h1 =subplot(2,2,1); 

% shade area of insignificance
x=[0 (nStim+2)*nElectrodes]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
shadedplot(x, y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on

% plot axis
for i = 1:9
    for j=1:nStim
        h3(j) =bar(j+(nStim+2)*(i-1), AzAll(i,j), 'FaceColor', barcol{j}, 'EdgeColor',barcol{j}*1.1,'LineWidth',.1);hold on
    end
%     if i == 9
%         plot(x, y1,'Color',[.25 .25 .25])
%         h1 =plot(4:9:76,mean(AzAll(:,1:4),2), '--b', 'LineWidth', 2);
%         h2 =plot(4:9:76,mean(AzAll(:,5:end),2), '--r', 'LineWidth', 2);
%     end
end

% axis and legend labels
set(gca, 'Xtick', 2:nStim+2:(nStim+2)*nElectrodes ,'xticklabel',{'EOG' 'Fp/aF' 'F' 'FC' 'C' 'CP' 'P' 'PO' 'O/I'},'FontSize', 8,'Ytick', 0:.25:1);
XMIN = 0;XMAX = (nStim+2)*nElectrodes-1;YMIN = .5;YMAX = 1;axis([XMIN XMAX YMIN YMAX]);ylabel('Az');
plot([0 (nStim+2)*nElectrodes-1],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1.3,'Color', [.25 .25 .25])

l2 = legend([h3(1:3)],{'BYD' 'GBU' 'PM'});
rect = [0.54, 0.745, .05, .1];set(l2, 'Position', rect);legend BOXOFF;
saveas(figure(4),'figure5','epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure 5- Free Viewing vs Constrained Viewing 
%%%% 
% note - include stimuli that are significant
homeDir = '/home/jason/Repository/AttentionModulation/';
run([homeDir 'Codes/eegInfo.m']) 
nModule = [1 2 3];

processIndx = 2 % choose manual or automatic processed file
processType = eeginfo.preprocess{processIndx}
iscFile =  [processType(1:end-1) '_allISC_v2.mat'];
saveDir = 'Data/ISCValues/'
saveDir = [homeDir saveDir];
load([saveDir iscFile]);

modality = [1 1 1 1 2 2 2 2 ];
attention = [1 2 1 2 1 2 1 2];
xAxis = [1 2 3 4 6 7 8 9 10];
narrative = [1 1 0 0 1 1 0 0];
subjgroup = [1 1 1 1 2 2 2 2];

% initialize variables
iscAll = []; iscSumAll = []; iscByCompAll =[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
iscCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];
iscSumMean =[];
compInt = 1:3 % Components of Interest
nComp = length(compInt) % measure first Ncomps

% plot/ figure settings
figure(5); clf;
barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
stimLabel = [1 1 2 2 1 1 2 2];
sp1= subplot(5,2,[1 3 5]);  hold on
stimOrder  = [1 3 2 4 5 7 6 8];
nStim = 8;
nSubjIndx=[1 1 1 1 2 2 2 2]

for iStim = 1:nStim
    isc_ = ISC{stimOrder(iStim)}(compInt,:);
    iscSum = sum(isc_,1);
    
    if iStim == 1; iscSum(end) =[];end
    iscSumMean(iStim) = mean(sum(isc_,1));
    iscSumAll = [iscSumAll; iscSum']; %con
     nSubj=size(iscSum,2);
    % concat labels
    subj = [1:nSubj]+(nSubjIndx(iStim)-1)*20 ;
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
    h(iStim)= bar(xAxis(iStim),mean(iscSum),0.5,'FaceColor', barcol{stimLabel(iStim)}, ...
        'EdgeColor',barcol{stimLabel(iStim)}*barColorScale,'LineWidth',.75,'LineStyle',line,'BarWidth',.75);
    plot(xAxis(iStim),iscSum,'.')
    yheight1(iStim) = mean(iscSum);
    n_{iStim} = num2str(length(iscSum));
end

nesting = [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
nesting = [0 1 0; 0 0 0 ; 0 0 0];

[P1,T1,STATS1,TERMS1] = anovan(iscSumAll,[subjLabelAll modalityIndxAll attentionIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1], ...
    'varname', {'subjects' 'visual condition' 'attention'}, 'display', 'on');


[P1,T2,STATS1,TERMS1] = anovan(iscSumAll(1:80),[subjLabelAll(1:80) attentionIndxAll(1:80)], ...
    'model', 'full', ...
    'nest', [], ...
    'random',[],...
    'varname', {'attention'}, 'display', 'off');

[P1,T3,STATS1,TERMS1] = anovan(iscSumAll(81:end),attentionIndxAll(81:end), ...
    'model', 'full', ...
    'nest', [], ...
    'random',[],...
    'varname', {'attention'}, 'display', 'off');

% initialize variables
iscAll = []; iscSumAll = []; iscByCompAll =[];
subjLabelAll= [];stimIndxAll=[];modalityIndxAll= [];
iscCompLabel = [];attentionIndxAll=[]; narrativeIndxAll=[];
iscSumMean =[];
compInt = 1:3 % Components of Interest
nComp = length(compInt) % measure first Ncomps 

attention = [1 2 1 2];
modality = [1 1 2 2];
sp1= subplot(5,2,[1 3 5]);  hold on
stimOrder  = [1 2 5 6]; %attend only
% stimOrder  = [3 4 7 8]; %count only
nStim = 4;
stimLabel=[1 2 1 2];
for iStim = 1:nStim
    isc_ = ISC{stimOrder(iStim)}(compInt,:);
    iscSum = sum(isc_,1);
    if iStim == 1; iscSum(end) =[];end
    iscSumMean(iStim) = mean(sum(isc_,1));
    iscSumAll = [iscSumAll; iscSum']; %concat iscs
    nSubj =  length(iscSum);
    nsubj_(iStim)=nSubj;
    % concat labels
    subj = [1:nSubj]+((subjgroup(iStim)-1)*size(iscSum,2)) ; 
    subjLabelAll = [subjLabelAll; subj'];

    stimIndx = repmat(stimLabel(iStim), nSubj,1); 
    stimIndxAll = [stimIndxAll; stimIndx];
    
    modalityIndx = repmat(modality(iStim), nSubj,1); 
    modalityIndxAll = [modalityIndxAll; modalityIndx];
    
end
meanFree=mean(iscSumAll(modalityIndxAll==1))
meanFixation=mean(iscSumAll(modalityIndxAll==2))
nesting = [0 1 0; 0 0 0; 0 0 0] 
[P1,T1,STATS1,TERMS1] = anovan(iscSumAll,[subjLabelAll modalityIndxAll stimIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1], ...
    'varname', {'subject' 'viewing condition' 'stimulus'}, 'display', 'on');

nesting = [0 1 0; 0 0 0; 0 0 0] 
[P2,T2,STATS2,TERMS2] = anovan(iscSumAll,[subjLabelAll modalityIndxAll stimIndxAll], ...
    'model', 'full', ...
    'nest', nesting, ...
    'random',[1], ...
    'varname', {'subject' 'viewing condition' 'stim'}, 'display', 'off');

clear p
stimIndx = [1 2 5 6; 3 4 7 8];
for iStim = 1:4
    if iStim==1;ISC{iStim}(:,end)=[];end
    attend = sum(ISC{stimIndx(1,iStim)}(compInt,:),1);
    disattend = sum(ISC{stimIndx(2,iStim)}(compInt,:),1);
    isc = [attend disattend]';
    attentionIndx= [ones(length(attend),1) ;zeros(length(disattend),1)];
    [h1(iStim),P{iStim}] = ttest(attend,disattend);
%     [P{iStim} t] = anovan(isc,attentionIndx,'display','off');
%     F(iStim) = t{2,6};
end

pAttend= .0002;
pDisattend= .0094;
lineStyle = {'-k','-k','-k','-k','-g','-g', '-g','-g' };
Y = [0 0 0 0 .17 .19 .17 .19];
p =mystack(cell2mat(P));
p(p > .05) = nan;
pBonf = bonf_holm(p);
pBonf(pBonf > .05) = nan
p = [pBonf;pAttend;pDisattend;1;1]
pFdr = mafdr(p,'BHFDR',1);
pFdr(pFdr > .05) = nan;
sigstarCustom({[1,2],[3,4], [6,7], [8,9], [1 8],[2 9], [3 6], [4 7]},p,[],lineStyle,Y);
set(sp1, 'Xtick', [2.5 7.5],'xticklabel', {'free viewing' 'constrained viewing'},'Ytick',[0:.05:.2]);
ylabel('ISC');box off
l = legend(h(1:4),{'BYD attend', 'BYD count', 'GBU attend', 'GBU count'})
rect = [0.5, 0.77, .02, .02]; set(l, 'Position', rect);legend BOXOFF
axis([0 10 0 0.2])
saveas(figure(5),'figure4','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 6- ORDER EFFECT COMPARISON (BYD,GBU,PM)
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
figure(6); clf; sp1= subplot(5,5,[1 2 3 6 7 8 11 12 13]); 

plotIndx  = [1 1 1 1 1 1 1 1 1 1 1 1]; 
subjgroup = [1 2 1 2 3 2 1 2 1 2 2 3];
stimOrder = [1 16 2 18 9 20 15 3 17 4 19 12];
attention = [1 2 1 2 1 2 1 2 1 2 1 2];

nStim = length(stimOrder); deleteNegativeValues =0;
iPlot = 1;
figure(6)
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
ylabel('ISC') ;axis([0 14 0 .215]);
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


saveas(figure(6),'figure6','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 7- ALPHA POWER ANALYSIS (BYD,GBU,PM,PMsc,Jpn)
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
f=figure(7);
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

% get Az 
attention =[1 2 9 10 11; 3 4 12 13 14];
stimOrder= [1 3 6 8 10; 2 4 7 9 11];
azIndx = [1 2 4 5 6];

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

subplot(2,9,[10:16])
% shade area of insignificance
x=[0 (nStim+2)]; y1= [sigPvalue sigPvalue] ; y2 = [0 0];
shadedplot([0 7], y1, y2,[ .7 .7 .7],[.25 .25 .25]); hold on



barcol = {[.3 .9 .9];[ 0 .6 .6]; [.5 .2 .8]; [.6 .4 .5]; [.7 .6 .5]};
for iStim =1:5
    hold on
    powerAttend = sum(power{attention(2,iStim)}(compInt,:),1);
    powerCount = sum(power{attention(1,iStim)}(compInt,:),1);
    label = [ones(length(powerAttend),1); zeros(length(powerCount),1)];
    azP(iStim) = rocarea([powerAttend powerCount], label);
    bar(azIndx(iStim), azP(iStim),0.5,'FaceColor', barcol{iStim}, ...
         'EdgeColor',barcol{iStim}*.02,'LineWidth',1,'BarWidth',.75)
end

plot([0 7],[sigPvalue sigPvalue],'LineStyle','-','LineWidth',1,'Color', [.25 .25 .25])
ylabel('Az');stims = {'BYD' 'GBU' 'PM' 'PMsc' 'Jpn'};
set(gca, 'Xtick', azIndx ,'xticklabel', stims, 'Ytick', [.5 .75 1], 'Yticklabel', {'0.5' '.75' '1'});
axis([0 7 .5 1])
saveas(figure(7), 'figure8', 'epsc')

% compare ROC
for iStim =1:5
    hold on
    iscAttend = sum(ISC{attention(1,iStim)}(compInt,:),1);
    iscCount = sum(ISC{attention(2,iStim)}(compInt,:),1);
    label = [ones(length(iscAttend),1); zeros(length(iscCount),1)];
    az_ISC(iStim) = rocarea([iscAttend iscCount], label);

    powerAttend = sum(power{attention(2,iStim)}(compInt,:),1);
    powerCount = sum(power{attention(1,iStim)}(compInt,:),1);
    label = [ones(length(powerAttend),1); zeros(length(powerCount),1)];
    az_Power(iStim) = rocarea([powerAttend powerCount], label);
    [pp{iStim} ,Ltheta{iStim},theta{iStim},S{iStim}] =comp2ROC(iscAttend,iscCount,powerAttend,powerCount)
end
