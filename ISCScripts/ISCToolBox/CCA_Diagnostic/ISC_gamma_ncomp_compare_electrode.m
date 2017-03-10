% Generate Figure 1
% BYD and GBU freeview attend vs disattend viewings 
% Generate comparison chart plus an ROC curve
% fileset{1 = Freeview attend and Freeview Disattend, 2 = Fixated attend 
%               and Fixated dissatend, 3 = Audio attend vs audio disattend}

clear all
close all
filedetails = '/home/jason/Experiments/Backwardcount/codes/bdf_reader/eegdetail.m';
run(filedetails)
eog = 0;


fileset = 3;
experiment = {    'freeview'    'fixated'    'audio'    };
  stims=Stims{fileset};stimsave=Stimsave{fileset};
    subjname=Subjname{fileset};eogchannel=Eogchannel{fileset};
    attend = Attend{fileset}; disAttend = Disattend{fileset};
load(['individual_iscs_by_electrodes_v0/' experiment{fileset} '.mat'], 'iscs');

Nstim = 3
nComp = 8;
nElectrode = 8
stim = [1:Nstim];
for iElectrode = 1:nElectrode
for iComp = 1:nComp
for iGamma = 1:size(iscs{1}.free{iElectrode}(:,:,:),3);
        sumisc1 = [];
        sumisc2 = [];
    for iStim = 2:Nstim
        isc_free = iscs{iStim}.free{iElectrode}(:,:,iGamma);
        isc_count = iscs{iStim}.count{iElectrode}(:,:,iGamma);
        Az(iElectrode,iGamma,iComp,iStim) = rocarea([sum(isc_free(1:iComp,:),1) sum(isc_count(1:iComp,:),1)], ...
            [ones(1,length(sum(isc_free(1:iComp,:),1))) zeros(1,length(sum(isc_count(1:iComp,:),1)))])
%         sumisc1 = [sumisc1 sum(iscs{iStim,1}.free(1:nComp,:),1)'];
%         sumisc2 = [sumisc2 sum(iscs{iStim,1}.count(1:nComp,:),1)'];
    end
    
%     sumisc = [sumisc1(:,stim) sumisc2(:,stim)];
%     sigtest1 = rm_anova2(sumisc(:),repmat((1:size(sumisc,1))',size(sumisc,2),1), ...
%         mystack(repmat([1:size(sumisc,2)/2 1:size(sumisc,2)/2],size(sumisc,1),1)), ...
%         mystack(repmat([1 2],(size(sumisc,2)/2)*size(sumisc,1),1)),{'movie','condition'});
%     [h p] = ttest(sumisc1(:,1),sumisc2(:,1));
%     [p h] = signtest(sumisc1(:,1),sumisc2(:,1));
%     Az;
end
end
end

 
% stims= {'pie Org' 'pie Scrm' 'Jpn'};

% for iStim =1:length(gamma)
% gammas{iStim} = num2str(gamma(iStim))
% end
for iStim = 1:Nstim
    figure(iStim)
    for i = 1:8
        
        subplot(8,1,i)
            bar(Az(:,:,i,iStim)')
            title(['Ncomp=' num2str(i)])
            axis([0 9 .5 1])
    end
        suptitle(stims{iStim})
    legend({'0' '.2' '.4' '.6' '.8'})

end

