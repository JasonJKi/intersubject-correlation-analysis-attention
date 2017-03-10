function data = harmonizeSampleRate(data,fsref,fs,N)
% harmonize sampling rate
% data = harmonizeSampleRate(data,fsref,fs,N)
    if fs~=fsref % Resample so that fs=fsref
        data2= resample(data,fsref,fs);
        if fs>fsref % Alter the trigger so that it has the correct number of values
            data2(:,N) = downsample(data(:,N),fs/fsref);
        else
            tmp = repmat(data(:,N),[1 fsref/fs])';
            data2(:,N) = tmp(:);
        end
        data=data2; clear data2;
    end
end
