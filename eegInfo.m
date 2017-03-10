%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess Information

eeginfo.moduleName = {'Freeview','Fixation','Audio','OrderReversed'};
eeginfo.eegProcDir = ['EEG/Processed/'];
eeginfo.eegUnprocDir = ['EEG/Unprocessed/'];
eeginfo.preprocess = {'Manual/' 'RPCA/' 'RPCA2/'}
eeginfo.stimName = {
    'BNG'
    'GBU'
    'PM'
    'PM scrm'
    'Jpn'
    };
eeginfo.stimIndx = [1 2 3 4 5];
eeginfo.fs = 256;
eeginfo.nChannels = 64;

eeginfo.stimDuration = [
    372
    400
    440
    444
    395
    ];

eeginfo
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment module information
module.eegStimName{1} = {
    'bng_1'
    'gbu_1'
    'bng_2'
    'gbu_2'
    };

module.outPutName{1} = { 
    'BYD free' 
    'GBU free'
    'BYD count' 
    'GBU count'
    };

module.subjRmv{1} = [5 19];
module.eogChannel{1} = 65:67;
module.attentionIndx{1} = [0 0 1 1];
module.stimIndx{1} = [1 2 1 2];
module.nSubj{1} = [23 22 22 22];

% Fixation Viewing Video
module.eegStimName{2} = {
    'bng_1'
    'gbu_1'
    'bng_2'
    'gbu_2'
    };

module.outPutName{2} = {
    'BNG fixate'
    'GBU fixate'
    'BNG fixate and count'
    'GBU fixate and count'
    };

module.subjRmv{2} = [];
module.eogChannel{2} = 65:70;
module.attentionIndx{2} = [0 0 1 1];
module.stimIndx{2} = [1 2 1 2];
module.nSubj{2} = [19 19 19 19];

% Free Viewing Audio
module.eegStimName{3} = {
    'pieOrg_1'
    'pieScrm_1'
    'jpn_1'
    'pieOrg_2'
    'pieScrm_2'
    'jpn_2'
    };

module.outPutName{3} = {    
    'PM free'
    'PM scrm. free'
    'Jpn free'
    'PM counting'
    'PM scrm. counting'
    'Jpn counting'};

module.subjRmv{3} = [];
module.eogChannel{3} = 65:70;
module.attentionIndx{3} = [0 0 0 1 1 1];
module.stimIndx{3} = [3 4 5 3 4 5];
module.nSubj{3} = [20 20 20 20 20 20];

% experiment module information
module.eegStimName{4} = {
    'bng_1'
    'bng_2'
    'gbu_1'
    'gbu_2'
    'pieOrg_1'
    'pieOrg_2'
    };

module.outPutName{4} = {
    'BYD count' 
    'BYD free' 
    'GBU count'
    'GBU free'
    'PM free'
    'PM count'
    };

module.subjRmv{4} = [];
module.eogChannel{4} = 65:70;
module.attentionIndx{4} = [0 0 1 1];
module.stimIndx{4} = [1 1 2 2 3 3];
module.nSubj{4} = [20 20 20 20];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/home/jason/AttentionModulation'))


