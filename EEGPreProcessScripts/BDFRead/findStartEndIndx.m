function [stindx endindx T M] = findStartEndIndx(trigger);
%  [stindx endindx T M] = findStartEndIndx(trigger);
% Find start and end indices for all trigger related events. 
% N is the trigger channel
% stindx = starting indices
% endindx = end indices
% T = last trigger recorded for the event, indicating duration of event.
% M = number of trigger events

    %Fix Errors in the trigger
    strg = 1
    j = 2; jj = 1;
    for i = 1:length(trigger)
        if i > 1 && trigger(i) ~= 0 && trigger(i-1) - trigger(i) < 0 ;
            sttrg(j) = trigger(i);
            if j > 2 && sttrg(j-1) ~= 255 && sttrg(j-1) ~= 69 && sttrg(j) - sttrg(j-1)  > 1;
                k = i;
                while trigger(k) == sttrg(j)
                    trigger(k) = 0;
                    k = k+1;
                end
            end
            j = j +1;
        end
    end
    
    % detect start (=0) and end times (coded as negative values)
    ttime = find(diff(trigger)>0)+1; % indices of triggers
    tvalue = trigger(ttime); % all trigger values
  
    
    
    % Case where tvalue(1) is not 69
    if tvalue(1) ~= 69
        trigger = [zeros(100, size(trigger,2)); trigger];
        ttime = find(diff(trigger)>0)+1; % indices of triggers
        tvalue = trigger(ttime); % all trigger values
    end
    
    tvalue = [tvalue; 0];
    for j=1:length(tvalue)-1
        
        % detect starts
        if tvalue(j)==69 && tvalue(j+1)==1
            tvalue(j) = 0;
        end
            
        % detect the end if
        if tvalue(j)>tvalue(j+1) && tvalue(j+1)~=1 ...% a jump down but not carryover
                || tvalue(j)==69 && tvalue(j+1)==69 ... % odd case of 69 sec long clip that is not the last one
                || tvalue(j)+1<tvalue(j+1) &&  tvalue(j+1)==69 % a jump up to by more than one to 69
            if tvalue(j-tvalue(j)) == 255
                tvalue(j) = tvalue(j)+256;
            end
            tvalue(j) = -tvalue(j);
        end
        j = j+1;
    end
    
    tvalue(find(tvalue == NaN)) = [];
    tvalue(end) = [];
    stind=find(tvalue==0);
    
    stindx = ttime(stind);
    
    % stindx = ttime(stind+1)-median(diff(ttime(stind*ones(1,100)+repmat((2:100+1),size(stind,1),1)),[],2)')';
    % % yuck!
    
    endindx1 = find(tvalue < 0);
    endindx = ttime(endindx1);
    
    T = abs(tvalue(find(tvalue < 0)));
    M = length(T);
end