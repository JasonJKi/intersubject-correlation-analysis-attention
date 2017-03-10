function timeParsedData = timeparseTriggerData(unParsedData, trigger, startTime, endTime,T, M)
%timeParsedData = timeparseTriggerData(unParsedData, trigger, startTime, endTime,T, M)
        
        x = unParsedData(startTime:endTime,:);
        trigger = trigger(startTime:endTime,:);
        
        % erase initial "69" of the trigger channel
        trg = trigger(:,end); 
        indx1=find(trg==1);
        trigger(1:indx1(1),end)=0; 
        trign = trigger(:,end); 
        plot(trign);
        
        % find all time indeces by seconds
        t=zeros(T,1);
        for n = 1:T
            tmp=find(rem(n,256)== trign);
            if n<256
                t(n)=tmp(1);                 %time markers
            elseif n==256
                t(n)=0;                      % fill in later;
            else
                tindx = tmp( find(diff(tmp)~=1)+ 1 );
                t(n)=tindx(end);
            end
        end
        
        % fix trigger at 256th frame for evengs longer then 256 seconds
        if length(t) > 256
            t(256)=round(mean([t(255) t(257)]));
        end
        
        % Perform time estimation for each trigger event
        A=[ones(T,1) [1:T]'];          % "mixing matrix"
        tmp=abs(pinv(A)*t);                  % estimated parameters
        t0=tmp(1);                           % estdispimated start time
        delta=tmp(2);                        % estimated sampling rate
        
        % estimated trigger event
        selrange = round(t0):round(t0)+T*round(delta);
        if min(selrange) == 0
            selrange = selrange+1;
        end

        timeParsedData = x(selrange,:);
end