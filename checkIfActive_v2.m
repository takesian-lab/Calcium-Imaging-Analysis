function [isActive, activity] = checkIfActive_v2(y, nBaselineFrames, STDlevel, AUC_level, plotFigure)

    baseline = y(1,1:nBaselineFrames);
    peak_threshold = nanmean(baseline) + STDlevel*std(baseline);
    trough_threshold = nanmean(baseline) - STDlevel*std(baseline);
    response = y(1,nBaselineFrames+1:end);
    
    if sum(baseline) == 0 %No activity in baseline
        isActive = 0;
        activity = 'none';
        return
    end

    %PEAK COMPUTATIONS
    peak_data = nan(1,4);
    [peak, peak_latency] = max(response);
    if peak >= peak_threshold %only store data if peak is above threshold
        [p1_latency] = find(response >= peak_threshold, 1, 'first');
        [p2_latency] = find(response(1, peak_latency:end) <= peak_threshold, 1, 'first') - 2;
        p1 = response(p1_latency);
        p2 = response(peak_latency + p2_latency);

        %AUC
        if isempty(p2_latency)
            p2_latency_temp = length(response);
        else
            p2_latency_temp = p2_latency + peak_latency;
        end
        peak_trace = response(1,p1_latency:p2_latency_temp);
        peak_trace(peak_trace < peak_threshold) = peak_threshold;
        peak_trace_no_nan = peak_trace(~isnan(peak_trace)); %trapz function does not work on nas
        aup = trapz(abs(peak_trace_no_nan - peak_threshold)); %Area under peak above threshold

        %Adjust for baseline
        peak_latency = peak_latency + nBaselineFrames;
        p1_latency = p1_latency + nBaselineFrames;
        p2_latency = p2_latency + peak_latency;

        %Width
        peak_width = p2_latency - p1_latency;
    else
        peak = nan;
        p1 = nan;
        p2 = nan;
        p1_latency = nan;
        p2_latency = nan;
        peak_latency = nan;
        peak_width = nan;
        aup = nan;
    end

    %Store
    if ~isempty(peak);          peak_data(1) = peak;            else;   peak = nan;         end
    if ~isempty(p1_latency);    peak_data(2) = p1_latency;      else;   p1_latency = nan;   end
    if ~isempty(peak_latency);  peak_data(3) = peak_latency;    else;   peak_latency = nan; end
    if ~isempty(peak_width);    peak_data(4) = peak_width;      else;   peak_width = nan;   end

    %TROUGH COMPUTATIONS
    trough_data = nan(1,4);
    [trough, trough_latency] = min(response);
    if trough <= trough_threshold %only store data if trough is below threshold
        [t1_latency] = find(response <= trough_threshold, 1, 'first');
        [t2_latency] = find(response(1, trough_latency:end) >= trough_threshold, 1, 'first') - 2;
        t1 = response(t1_latency);
        t2 = response(trough_latency + t2_latency);

        %AUC
        if isempty(t2_latency)
            t2_latency_temp = length(response);
        else
            t2_latency_temp = t2_latency + trough_latency;
        end
        trough_trace = response(1,t1_latency:t2_latency_temp);
        trough_trace(trough_trace > trough_threshold) = trough_threshold;
        trough_trace_no_nan = trough_trace(~isnan(trough_trace));
        aat = trapz(abs(trough_trace_no_nan - trough_threshold)); %Area above trough and below threshold

        %Adjust for baseline
        trough_latency = trough_latency + nBaselineFrames;
        t1_latency = t1_latency + nBaselineFrames;
        t2_latency = t2_latency + trough_latency;

        %Width
        trough_width = t2_latency - t1_latency;
    else
        trough = nan;
        t1 = nan;
        t2 = nan;
        t1_latency = nan;
        t2_latency = nan;
        trough_latency = nan;
        trough_width = nan;
        aat = nan;
    end

    %Store
    if ~isempty(trough);            trough_data(1) = trough;            else;   trough = nan;           end
    if ~isempty(t1_latency);        trough_data(2) = t1_latency;        else;   t1_latency = nan;       end
    if ~isempty(trough_latency);    trough_data(3) = trough_latency;    else;   trough_latency = nan;   end
    if ~isempty(trough_width);      trough_data(4) = trough_width;      else;   trough_width = nan;     end

    %Auto-determined activity (inhibited/sustained/activated)
    if ~isnan(aup)&& aup >= AUC_level; aup_pass = true; else; aup_pass = false; end
    if ~isnan(aat)&& aat >= AUC_level; aat_pass = true; else; aat_pass = false; end

    activity = 'undetermined'; %If it somehow makes it through the conditions without being classified

    if isnan(peak) && isnan(trough)
        activity = 'none';
    elseif ~aat_pass && ~aup_pass
        activity = 'none';
    elseif isnan(peak) && ~isnan(trough) && aat_pass
         activity = 'inhibited';
    elseif ~isnan(peak) && isnan(trough) && aup_pass
        if peak_latency > 40 || isempty(p2_latency)
            activity = 'sustained';
        else
            activity = 'activated';
        end
    elseif ~isnan(peak) && ~isnan(trough)
        if (trough_latency < peak_latency) && aat_pass
            activity = 'inhibited';
        elseif (peak_latency < trough_latency) && aat_pass && ~aup_pass
            activity = 'inhibited';
        elseif aup_pass
            if peak_latency > 40 || isempty(p2_latency)
                activity = 'sustained';
            else
                activity = 'activated';
            end
        else
            activity = 'none';
        end
    else
        activity = 'none';
    end

    if ~strcmp(activity, 'undetermined') && ~strcmp(activity, 'none')
        isActive = 1;
    else
        isActive = 0;	
    end
    
    %Plot
    if plotFigure
        figure; hold on
        plot(y)
        hline(nanmean(baseline), 'k')
        hline(peak_threshold, 'r')
        hline(trough_threshold, 'c')
        scatter(peak_latency, peak, 'o', 'r')
        scatter(p1_latency, p1, 'o', 'r')
        scatter(p2_latency, p2, 'o', 'r')
        scatter(trough_latency, trough, 'o', 'c')
        scatter(t1_latency, t1, 'o', 'c')
        scatter(t2_latency, t2, 'o', 'c')
        vline(nBaselineFrames, 'k')
        xlabel('Frames')
        %ylabel(units)
        if strcmp(activity, 'activated') || strcmp(activity, 'sustained')
            title([activity ' -  AUC: ' num2str(aup)])
            plot(p1_latency:(nBaselineFrames + p2_latency_temp), peak_trace, 'g')
        elseif strcmp(activity, 'inhibited')
            title([activity ' - AUC: ' num2str(aat)])
            plot(t1_latency:(nBaselineFrames + t2_latency_temp), trough_trace, 'g')
        else
            title(activity)
        end
    end

    end