function [eegclean samplesRemoved] = artifactRejection(eeg,eogchannels,nChannels, fsref)
% [eegclean samplesRemoved] = preprocessingThreshRemoval(eeg,eogchannels,nChannels, fsref)
% 
% preprocess eeg and label eeg quality
% steps
% 1) EOG regression 
% 2) Artifact rejection by thresh holding using standard deviation across
% channels.
% 3)  Artifact rejection by thresh holding using standard deviation across
% time samples.
Sremove = 3; % number standardation removal
Nremove = 4; % iteration of removal

% Electrodes
EOG = eeg(:,eogchannels);

% EEG Electrodes = 1-64
EEG = eeg(:,1:nChannels);


% remove eye movement artefact as best as we can
eegclean = EEG - EOG * (EOG\EEG);
% figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1); imagesc(EEG(fsref*10+(1:fsref*10),:)'); caxis([-100 100]); title('bad channel removed')
subplot(2,3,4); imagesc(eegclean(fsref*10+(1:fsref*10),:)'); caxis([-100 100]);title('eye movement regressed')

% artifact rejection channel
mask = zeros(size(eegclean));
se = ones(round(fsref*0.025),1);
for k=1:Nremove;
    thresh = Sremove*std(eegclean);
    mask = mask + filtfilt(se,1,double(abs(eegclean)>repmat(thresh,[size(eegclean, 1) 1])));
    eegclean(mask>0)=0;
end
maskc =mask'>0;

k = repmat(thresh',[1 size(eegclean, 2)]);
subplot(2,3,2); imagesc(eegclean'); caxis([-100 100]);title('artifact removed')
subplot(2,3,5); imagesc(maskc); title('time artifacts')

% artifact rejection time sample
mask = zeros(size(eegclean));
se = ones(round(fsref*0.01),1);
for k=1:2
    thresh = Sremove*std(eegclean');
    mask = mask + filtfilt(se,1,double(abs(eegclean)>repmat(thresh',[1 size(eegclean, 2)])));
    eegclean(mask>0)=0;
end
maskt = mask'>0;

subplot(2,3,3); imagesc(eegclean'); caxis([-100 100]);title('artifact removed')
subplot(2,3,6); imagesc(maskt); title('channel artifacts')

samplesRemoved = sum(maskc(:)) + sum(maskt(:))
suptitle(['num samples removed = ' num2str(samplesRemoved)])