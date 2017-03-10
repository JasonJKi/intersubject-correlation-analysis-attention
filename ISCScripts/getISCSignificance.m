% compute the ISC significance.
clear all;
homeDir = '/home/jason/AttentionModulation/';
run([homeDir 'eegInfo.m'])
nModule = [1 2 3];
processtype = eeginfo.preprocess{2}
currentDir = 'ISCSignificance/'
compIndx = [1 1 2]; nComp = 10;

i = 1;
for iModule = nModule
    nStim = length(module.eegStimName{iModule});
    eogIndx = module.eogChannel{iModule};
    moduletype = eeginfo.moduleName{iModule}
    for iStim = 1:nStim;
        stim = module.eegStimName{iModule}{iStim}
        eeginDir = [homeDir eeginfo.eegProcDir processtype moduletype '/' stim '.mat'];
        load(eeginDir)
        eeg(:,:,module.subjRmv{iModule}) = [];
        [t nChannels nSubj] = size(eeg);
        for iPermutation = 1:100;
            
            % method 1 - online
            T = t;
            if rem(t,2)==0;T = T-1;
            end
            eegRandomizedAll_ = zeros(T, nChannels, nSubj);
            for iSubj = 1:nSubj
                disp(iSubj)
                % Fourier transform of the original dataset
                eegFFT = fft(eeg(1:T,:,iSubj));
                
                % Create the random phases for all the time series
                ph_interv1 = repmat(exp( 2*pi*1i*rand([lengthEEG 1])),1,nChannels);
                ph_interv2 = conj( flipud( ph_interv1));
                
                % Randomize all the time series simultaneously
                fft_recblk_surr(2:(T-1)/2+1,:) = eegFFT(2:(T-1)/2+1,:).*ph_interv1;
                fft_recblk_surr((T-1)/2+2:T,:) = eegFFT((T-1)/2+2:T,:).*ph_interv2;
                
                % Inverse transform
                eegRandomizedAll_(:,:,iSubj)= real(ifft(fft_recblk_surr));
            end
            [Rxy Rpool] = generate_cov(eegRandomizedAll_,0);
            
            
            % method 2 - Samantha
            T = t;
            if mod(T,2)~=0;T = T-1;end
            for iSubj = 1:nSubj
                disp(iSubj)
                % Fourier transform of the original dataset
                eegFFT = fft(eeg(1:T,:,iSubj)); 
                %compute frequency amplitude
                amplitude = abs(eegFFT(1:T/2+1,:)); 
                %compute phase
                phase = angle(eegFFT(1:T/2+1,:)); 
                % Generate random phase
                phaseRandomized = (-pi) + (2*pi)*rand(T/2-1,1);
                
                eegRandomized(2:T/2,:) = amplitude(2:T/2,:).*exp(1i*(phase(2:T/2,:)+repmat(phaseRandomized,1,nChannels))); % randomize the phase
                eegRandomized = [eegFFT(1,:); eegRandomized(2:T/2,:); eegFFT(T/2+1,:); conj(eegRandomized(T/2:-1:2,:))]; % rearrange  freq
                eegRandomizedAll(:,:,iSubj) = ifft(eegRandomized); % reverse back to time domain
            end
            disp(iPermutation)
            [Rxy Rpool] = generate_cov(eegRandomizedAll,1);
            [w a] = correlated_componenChannels(Rxy,Rpool,.5,64, 0);

            for iSubj = 1:nSubj
            isc(:,iSubj) = concat_matrix_ISC(eegRandomizedAll,eegRandomizedAll,w,iSubj,10);
            end
            
            iscChance(:,iPermutation) = mean(isc,2)
            
        end
        iscChanceAll(:,:,i) = iscChance;
        i = i + 1;
    end
end

