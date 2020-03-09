function [peakLocationfinal, peakHeight, peakRatio, troughLocationfinal, SNR, RMSE] = peakLocation_prominence_2019_update(movingmeanTrace, imagingspeed, sensor)
%Inputs
%movingmeanTrace: An optical trace that has been lowpass filtered
%imagingspeed: Speed in Hz of video acquisition. Temporally binned is 5Hz
%otherwise it's 20Hz
%sensor: GCaMP6f or GCaMP6s. Used to determine the criteria for identifying
%transients
%multipeak: If it's 1 it will identify all troughlocations within a given
%short bought of high DF/F activity, if it's 0 it will only identify the
%initial trough before DF/F activty began.
traceout = double(movingmeanTrace(:, 1));


%Identifying calcium transients presents an interesting problem. We need to
%identify the baseline noise of the signal, but there are far too many
%neurons being analyzed for this to be a problem that can be solved by
%brute force. Instead, we need to identify a way of automatically
%thresholding our data to identify calcium transients. The standard
%deviation is a poor choice because it is heavily influenced by the
%presence of calcium transients. Instead I use the median absolute
%deviation as it is a more robust measurement of the baseline noise of
%calcium activity. 
try
    medabsdev = mad(traceout, 1);
    baseMean = median(traceout);
    threshold = medabsdev*4;
    
    %finds all indexes of the trace above the threshold
    rind = find(traceout > (baseMean+threshold));
    %If there are no points above the threshold then there are no
    %transients to be identified in the trace
    if isempty(rind) == 1 || length(rind) == 1
        start = NaN(1);
        peakHeight = NaN(1);
        peakRatio = NaN(1);
        troughLocation = NaN(1);
        troughLocationfinal = NaN(1);
    else
        %Finds the start and end points for each short chunk of the trace that
        %is above threshold
        drind = diff(rind);
        changes = find(drind ~= 1);
        
        start = rind([1 changes(1:end)'+1]);
        ends = rind([changes' size(rind, 1)]);
        peakLocationfinal = [];
        troughLocationfinal = [];
        peakHeight = [];
        %Determines whether these short periods of time are longer than a set
        %threshold, and also exludes transients that occur at the beginning of
        %the video since I can't accurately assess when the trough occurs
        for fn = 1:size(start, 1)
            refind = start(fn):ends(fn);
            switch sensor
                case 'GCaMP6s'
                    minLength = imagingspeed;
                    risetime = ceil(imagingspeed*0.2);
                case 'GCaMP6f'
                    minLength = ceil(imagingspeed*0.3);
                    risetime = ceil(imagingspeed*0.2);
            end
            %Finds the first initial peak, and then sets a higher threshold for
            %any following peaks
            if length(refind) > minLength && start(fn) ~= 1
                shorttrace = traceout(start(fn):ends(fn));
                
                [clumppeaks, shortlocations] = findpeaks(shorttrace);
                touse = [];
                lastaccepted = [];
                %If there's more than one identified peak within the region
                %of fluorescence activity that was above baseline we need
                %to identify legitimate peaks. To do this we use a couple
                %of rules: 
                %1: the peak must be 2.5 median absolute deviations above
                %the next identified peak
                %2: The peak must be 5 median absolute deviations above the
                %previously identified peak
                %3: The maximum peak is always accepted in this algorithm
                
                if length(clumppeaks) > 1
                    for ss = 1:length(clumppeaks)
                        if ss == 1
                            if clumppeaks(ss) > clumppeaks(ss+1)+2.5*medabsdev
                                touse(ss) = 1;
                            else
                                touse(ss) = 0;
                            end
                            lastaccepted = 1;
                        elseif ss == length(clumppeaks)
                            if clumppeaks(ss) > clumppeaks(lastaccepted)+5*medabsdev
                                touse(ss) = 1;
                          
                            else
                                touse(ss) = 0;
                            end
                        else
                            if clumppeaks(ss) > clumppeaks(ss+1)+2.5*medabsdev && clumppeaks(ss) > clumppeaks(lastaccepted)+5*medabsdev
                                touse(ss) = 1;
                                lastaccepted = ss;
                            else
                                touse(ss) = 0;
                            end
                        end
                        if clumppeaks(ss) == max(clumppeaks)
                            touse(ss) = 1;
                        end
                    end
                    shortref = shortlocations(logical(touse));
                    allpeak = clumppeaks(logical(touse));
                else
                    shortref = shortlocations;
                    allpeak = clumppeaks;
                end
                
                touse2 = [];
                cutoffs = allpeak*exp(-1);
                %Now we do a quick check to make sure that the calcium
                %transient is not decaying faster than we would predict
                %with an exponential decay. This account for sudden shifts
                %that might occur due to a registration error. 
                for ss = 1:length(shortref)
                    if refind(shortref(ss))+ 3 <= length(traceout)
                    if traceout(refind(shortref(ss))+3) > cutoffs(ss)
                        touse2(ss) = 1;
                    else
                        touse2(ss) = 0;
                    end
                    else
                        touse2(ss) = 1;
                    end
                end
                shortref = shortref(logical(touse2));
                allpeak = allpeak(logical(touse2));
                
                %We have identified the peaks. Now we want to generate a
                %list of the troughs. For the first peak the trough is the
                %first frame that surpasses threshold. For all other peaks
                %it is the local minima between peaks. 
                if ~isempty(shortref)
                    traceLocations = refind(shortref);
                    peakLocationfinal = cat(1, peakLocationfinal, traceLocations');
                    peakHeight = cat(1, peakHeight, allpeak);
                    troughlocations = [];
                    for ss = 1:length(shortref)
                        if ss == 1
                            troughlocations(ss) = refind(1);
                        else
                            [minval minloc] = min(shorttrace(shortref(ss-1):shortref(ss)));
                            troughlocations(ss) = refind(shortref(ss-1)+minloc-1);
                        end
                    end
                    troughLocationfinal = cat(1, troughLocationfinal, troughlocations');
                    
                    
                end
            else
            end
        end
        
        
    end
    %Calculate the ratio of the peak size and the SNR of the transients
    peakRatio = peakHeight./min(peakHeight);
    SNR = (mean(peakHeight)-baseMean)/standardd;
catch
    disp('Coulndnt fit')
    peakLocationfinal = NaN(1);
    peakHeight = NaN(1);
    peakRatio = NaN(1);
    troughLocationfinal =NaN(1);
    SNR = NaN(1);
    RMSE = NaN(1);
end

end

