% Read in raw bdf and create matlab arrays.
clear all;
fileDir='/home/jason/Repository/AttentionModulation/Data/EEG/';

% Indicate desired save path for the unprocessed EEG Mat files.
exptype = 'OrderReversed2';
savepath1 = [fileDir 'Unprocessed/' exptype '/'];
disp(['File will be stored in         ' savepath '        press any key to continue'])
mkdir(savepath1)
% Biosemi - 64 Channels
Nchannel = 64;
% Desired Sampling Rate
fsref = 256;

% Lowpass/drift and 60Hz noise remove filter
[hpnum,hpdenom]=butter(4,0.3/fsref*2,'high'); % drift removal
[notchnum,notchdenom]=butter(4,[59 61]/fsref*2,'stop'); % 60Hz line noise
a = poly([roots(hpdenom);roots(notchdenom)]);
b = conv(hpnum,notchnum);

freqz(a,b)
filobj = [a; b];
savemethod = 'subject';

stimNames= {'bng' 'gbu' 'pieOrg'};
endTrigger= [
    373
    401
    441
    445
    396
    ];

stimDuration=[374.2800 389.0050 427.9253]
fsRatio=[1 29.97/29 29.97/29 29.97/29];
bdfPath = '/home/jason/Repository/AttentionModulation/Data/BDF/OrderReversed/'
for ii = 17:20;
%   BDFtoTriggerEvents(bdffile, savepath1, stimnames, duration, fsref, Nchannel, savemethod)
%   parse eeg trigger events by sample and movie title
%   savemethod = 'subject' or 'stim'
    subject = ['Subject' num2str(ii) '.bdf']
    BDFtoTriggerEvents(bdfPath, subject ,savepath1, stimNames, endTrigger, filobj, fsref,fsRatio,stimDuration);
end

