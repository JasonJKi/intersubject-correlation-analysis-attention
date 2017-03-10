% Generate Figure 1
% BYD and GBU freeview attend vs disattend viewings 
% Generate comparison chart plus an ROC curve
% fileset{1 = Freeview attend and Freeview Disattend, 2 = Fixated attend 
%               and Fixated dissatend, 3 = Audio attend vs audio disattend}

clear all
filedetails = '/home/jason/Experiments/Backwardcount/codes/bdf_reader/eegdetail.m';
run(filedetails)
eog = 0;


fileset =1;
experiment = {    'freeview'    'fixated'    'audio'    };
load(['individual_iscs_v2/' experiment{fileset} '.mat'], 'iscs');

Nstim = length(Stims{fileset})/2
nComp = 10;
stim = [1:Nstim]
for iComp = 1:nComp
for iGamma = 1:size(iscs{1}.free(:,:,:),3);
        sumisc1 = [];
        sumisc2 = [];
    for iStim = 1:Nstim
        isc_free = sum(iscs{iStim}.free(1:iComp,:,iGamma),1);
        isc_count =sum(iscs{iStim}.count(1:iComp,:,iGamma),1);
        Az(iComp,iGamma,iStim) = rocarea([isc_free isc_count], ...
            [ones(1,length(isc_free)) zeros(1,length(isc_count))])
         sumisc1 = [sumisc1 sum(isc_free,1)'];
         sumisc2 = [sumisc2 sum(isc_count,1)'];
    end
    
     sumisc = [sumisc1(:,stim) sumisc2(:,stim)];
     sigtest1{iComp,iGamma} = rm_anova2(sumisc(:),repmat((1:size(sumisc,1))',size(sumisc,2),1), ...
         mystack(repmat([1:size(sumisc,2)/2 1:size(sumisc,2)/2],size(sumisc,1),1)), ...
         mystack(repmat([1 2],(size(sumisc,2)/2)*size(sumisc,1),1)),{'movie','condition'});
     [h(iComp,iGamma) p(iComp,iGamma)] = ttest(sumisc1(:,1),sumisc2(:,1));
     [p(iComp,iGamma) h(iComp,iGamma)] = signtest(sumisc1(:,1),sumisc2(:,1));
     Az;
    end
    
end


 
% stims= {'pie Org' 'pie Scrm' 'Jpn'};
gamma = .1:.1:.9
for iStim =1:length(gamma)
gammas{iStim} = num2str(gamma(iStim));
end
figure(fileset)
for iStim = 1:Nstim
   subplot(Nstim,1,iStim)
    bar(Az(:,:,iStim),'grouped')
    ylabel(Aztitle{fileset}{iStim})
    xlabel(['# of comp' iComp])
    axis([0 11 .5 1 ])
end

    suptitle(experiment{fileset})
    legend(gammas)