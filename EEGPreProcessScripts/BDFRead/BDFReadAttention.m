% Read in raw bdf and create matlab arrays.

clear all;

% Indicate desired save path for the unprocessed EEG Mat files.
exptype = 'fixated';
mkdir('/home/jason/Experiments/Backwardcount/data/unprocessed/', exptype)
savepath1 = ['/home/jason/Experiments/Backwardcount/data/unprocessed/' exptype '/'];
disp(['File will be stored in         ' savepath '        press any key to continue'])

% Choose file set (Freeviewing = 1, Fixation = 2, Audio Only = 3)

% Open file and stim list
run EEGDetail.m;
fileset = 3
filename = filenames{fileset};
filepath = filepaths{fileset};% EEG specification
n = length(filename);
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

disp(filename)
disp(filepath)
savemethod = 'subject';

for ii = 19:20;
%   BDFtoTriggerEvents(bdffile, savepath1, stimnames, duration, fsref, Nchannel, savemethod)
%   parse eeg trigger events by sample and movie title
%   savemethod = 'subject' or 'stim'
    BDFtoTriggerEvents(filepath, filename{ii}, savepath1, stim.name, stim.duration, filobj, fsref);

end
