function[spike_peaks_filt,raw_peaks,spike_peaks,section,closest_to_start] = find_spikes(ASAPtime,cutoff,ratio,ff,fF,smooth_wave,data,ratio_thresh,updown,pre_spike,post_spike,fs,tr,varargin);

% After adjusting the peaks to shift so that the ratio peaks map to the
% peaks in the raw data, it is possible that 2 peaks get moved to the same
% location if two are detected within 2ms of each other. To prevent this, the duplicate peak(s) is discarded.

% Note that it chooses which spike position to keep based on the size of
% the largest one in the data vector.

% Currently spike_peaks_filt isn't lining up to the raw peaks correctly
% because we select both based on largest response
% in a window, and sometimes after filtering
% different ones have a different largest response.

[~,closest_to_start] = min(abs(ASAPtime-cutoff));
                temp = find(ratio(:,tr)>ratio_thresh);
                Above_thresh = temp(temp>=closest_to_start); % Inds where ratio is above thresh, minus cutoff
                section = {};
                spike_peaks = [];
                spikes = [];
                raw_spikes = [];
                ratio_range = [];
                spike_peaks_filt = [];
                raw_peaks = [];

                if isempty(Above_thresh)==0 % If there are spikes detected
                    
                    diff_above_thresh = diff(Above_thresh); % The difference between the inds where we are above thresh.
                    % If diff is 1, then it means 2 time points are part of
                    % the same event.
                    
                    logical_above_thresh = diff_above_thresh==1; % Clusters
                    % of ones correspond to a single spiking event, noted
                    % by the corresponding inds in Above_thresh. Note that
                    % the ind after the last 1 of each cluster is also
                    % included in the previous cluster, because of the way
                    % differencing works.    
                    
                    % The first 0 after a string of ones in
                    % logical_above_thresh is the last point of a long
                    % event above threshold. Any zeroes after that are
                    % single point events where an event only crossed
                    % threshold for a single data point.
                    
                    new_temp = [];
                    for single_inds = 2:length(logical_above_thresh)
                       
                        if logical_above_thresh(single_inds)==0 && logical_above_thresh(single_inds-1)==0
                            new_temp(single_inds) = 1; % Marks all of the locations of single spikes
                        else new_temp(single_inds) = 0;
                        end
                        
                    end

                    single_inds = find(new_temp>0);
                    
                    ratio_spike_peak_inds = Above_thresh(diff_above_thresh>1); % Using the point where ratio first crosses threshold to count as spike onset.
                    % Not great because it will vary dramatically depending
                    % on how the threshold is set.
                    
                    find_edges = evans_conv([-1,1],logical_above_thresh,1,0); % Gives a 0.5 when a 1 starts, and a -0.5 when it hits a 0.
                    start_inds_above_thresh = find(find_edges==0.5); 
                    
                    start_inds_above_thresh = [single_inds';start_inds_above_thresh'];
                    start_inds_above_thresh = sort(start_inds_above_thresh);
                    
                    for k1 = 1:length(start_inds_above_thresh) % 1 loop for each start detected

                        on_deck = start_inds_above_thresh(k1);
                        counter = 0;

                            while logical_above_thresh(on_deck+counter)==1
                                counter = counter+1;

                                if on_deck+counter>length(logical_above_thresh);
                                    break;
                                end % get out of the while loop. logical
                                % above_thresh ends on a 1.
                            end
                            end_inds_above_thresh(k1) = on_deck+counter; % Segment End Indices 
                            section{k1} = [Above_thresh(start_inds_above_thresh(k1):end_inds_above_thresh(k1))];
                            [~,spike_peaks(k1)] = max(ratio(section{k1},tr)); % Find ind of max within spiking event
                            spike_peaks(k1) = spike_peaks(k1)+Above_thresh(start_inds_above_thresh(k1))-1; % Adjust ind
                            % to fit into original vector..
                    end
                    if isempty(k1)==1
                        k1 = 0;
                    end

                % Add in the single point spikes
                spike_peaks = sort(spike_peaks);
                % Test figure to make sure it is grabbing appropriate peaks
                %figure;plot(ratio(:,tr));hold on;plot(spike_peaks,ratio(spike_peaks,tr),'r*');hold off                                       

                spike_peaks = spike_peaks(spike_peaks/fs>cutoff); % Chop beginning
                spike_peaks = spike_peaks(spike_peaks/fs<(length(ff)/fs)-cutoff); % Chop end
                
                      on_deck = [];
                      raw_locs = [];
                      raw_pks = [];
                      raw_peaks = [];
                      peaks = [];
                      on_deck_filt = [];
                      spike_peaks_filt = [];
                      pks = [];
                      locs = [];
                      spike_count_raw = 0;
                      spike_count_filt = 0;
                      search_width = fs/10; % Smaller = search a wider area for peak. This currently corresponds to a search area of 5 data points +/- from the peak at 1000hz
                      
                     for j = 1:length(spike_peaks) % 1 loop for each spike. Peak in trace isn't always the same
                         % as the peak in ratio, sometimes off by a few
                         % sample points, so we have to shift it.
                         
                         start_ind = spike_peaks(j)-round(fs/search_width);
                         if start_ind<1;start_ind = 1;end
                         
                         end_ind = spike_peaks(j)+round(fs/search_width);
                         if end_ind>length(data(:,tr))
                             end_ind = length(data(:,tr));
                         end
                         
                        on_deck = data(start_ind:end_ind,tr);
                        on_deck_filt = smooth_wave(start_ind:end_ind,tr);
                        
                        if updown == 1 % Depends on if the indicator goes up or down
                            [raw_pks(j), raw_locs(j)] = max(on_deck);
                            [raw_pks_filt(j), raw_locs_filt(j)] = max(on_deck_filt);                            
                        else
                            [raw_pks(j), raw_locs(j)] = min(on_deck);
                            [raw_pks_filt(j), raw_locs_filt(j)] = min(on_deck_filt);  
                        end                                              
                        
                        raw_peaks_temp = spike_peaks(j)-(round(fs/search_width)+1)+raw_locs(j);
                        filt_peaks_temp = spike_peaks(j)-(round(fs/search_width)+1)+raw_locs_filt(j);

                            spike_count_raw = spike_count_raw+1;
                            raw_peaks(spike_count_raw) = raw_peaks_temp; % Find biggest change in the spikes within the window                                             
                            spike_peaks_filt(spike_count_raw) = filt_peaks_temp;
                            
                     end % end for j loop

                     % If the raw peak is closest to the edge
                        raw_inds_hold = raw_peaks;
                        spike_peaks_filt_hold = spike_peaks_filt;
                        
                        raw_spike_peaks_temp = raw_peaks - pre_spike; raw_spike_peaks_temp = raw_spike_peaks_temp>0;
                        raw_peaks = raw_peaks(raw_spike_peaks_temp);
                        
                        raw_inds_hold = raw_inds_hold(raw_spike_peaks_temp);
                        spike_peaks_filt_hold = spike_peaks_filt_hold(raw_spike_peaks_temp);
                        
                        section = section(raw_spike_peaks_temp);
                        [raw_peaks,IA,IC] = unique(raw_peaks);
                        section = section(IA);
                        
                    % Now search for the corresponding raw peak in the ratio and filtered vectors.
                    
                    spike_peaks_temp = [];
                    filt_peaks_temp = [];
                    
                    for duplicates = 1:length(raw_peaks)
                       
                        [~,I] = min(abs(spike_peaks - raw_peaks(duplicates)));
                        spike_peaks_temp(duplicates) = spike_peaks(I);
                        
                        [~,I] = min(abs(spike_peaks_filt - raw_peaks(duplicates)));
                        filt_peaks_temp(duplicates) = spike_peaks_filt(I);
                        
                    end                   
                    spike_peaks = spike_peaks_temp;
                    spike_peaks_filt = filt_peaks_temp;
                                        
                    %% Now that we have the raw, ratio and filtered vector
                    % peaks, take out spikes near cutoff that would extend over. Do
                    % this by taking the ind closest to the edge.
                    temp = [spike_peaks(1) raw_peaks(1)];
                    [~,I] = min(temp);

                    spike_peaks_temp = []; % Reset the variable
                    if I==1 % If the ratio peak is closest to the edge
                        spike_peaks_temp = spike_peaks - pre_spike; spike_peaks_temp = spike_peaks_temp>0;
                    end
                    if I==2 % If the raw peak is closest to the edge
                        spike_peaks_temp = raw_peaks - pre_spike; spike_peaks_temp = spike_peaks_temp>0;
                    end

                    spike_peaks = spike_peaks(spike_peaks_temp);
                    raw_peaks = raw_peaks(spike_peaks_temp);
                    section = section(spike_peaks_temp);
                     
                end % end for if we found spikes                
                
if sum(strcmp(varargin{1},'min_interval'))>0 % If the user input a minimum interval
    % Remove events that are within min_interval of eachother and are
    % smaller than the close event
    
    temp = varargin{1,1};
    temp2 = find(strcmp(varargin{1},'min_interval')==1)+1;
    min_interval = str2double(string(temp(temp2)));    
    
    ratio_peaks_seconds = (spike_peaks/fs); % Convert ratio peaks into seconds
    ratio_peaks_diff = diff(ratio_peaks_seconds); % Seconds between events
    peaks_within_interval = ratio_peaks_diff<=min_interval; % Logical vector, a 1 means that event and the one in front are within min_interval of eachother
    find_segments = evans_conv([-1,1],peaks_within_interval,1,0); % Gives a 0.5 when a 1 starts, and a -0.5 when it hits a 0.
    
    % Note that segments correspond to the gap between events, so there
    % will be one less than the number of events
    %% Convert it to correspond to actual indices in the events:

    find_segments_new = [];
    combine_switch = 0;

for i = 1:length(find_segments)
   
    if find_segments(i)==0 & combine_switch==0 % If we find_segments_newren't combining events
        find_segments_new(i) = 0;
    end
    if find_segments(i)==0 & combine_switch==1
        find_segments_new(i) = 1;
    end
    if find_segments(i)==0.5
        find_segments_new(i) = 1;
        combine_switch = 1; % Stfind_segments_newrt combining sections
    end
    if find_segments(i)==-0.5
        find_segments_new(i) = 1;
        combine_switch = 0; % Stop combining sections
    end
        
end

if find_segments_new(end)==0 & find_segments(end)==0
    find_segments_new(length(find_segments_new)+1) = 0; % Not combining with lfind_segments_newst event
end
if find_segments_new(end)==0 & find_segments(end)==0.5 % If find_segments_new combine event stfind_segments_newrted find_segments_newt the end
    find_segments_new(length(find_segments_new)+1) = 1; % Combining with lfind_segments_newst event
end


if find_segments_new(end)==1 & find_segments(end)==-0.5
    find_segments_new(length(find_segments_new)+1) = 0; % Not combining with lfind_segments_newst event
end
if find_segments_new(end)==1 & find_segments(end)==0.5 % If find_segments_new combine event stfind_segments_newrted find_segments_newt the end
    find_segments_new(length(find_segments_new)+1) = 1; % Combining with lfind_segments_newst event
end


%% END convert it to correspond to actual event indices
    
    ratio_inds = [];
    
    combine_switch = 0; % Changes to a 1 when a 0.5 is detected (start of combine). Once it becomes
    % a -0.5 (end of combine these sections) it changes back to 0, which
    % means add these sections individually.
    new_events = 0; % This will count the new events, as we will be combining them so there will be fewer than original
    
    section_hold = [];
    spike_peaks_filt_hold = [];
    raw_peaks_hold = [];
    spike_peaks_hold = [];
    spike_peaks_temp = [];
    
    for join = 1:length(find_segments_new) % Combine sections that are too close. 0.5 means combine this section
        % and the rest until you hit a -0.5. 
        
        if find_segments_new(join)==1 % If we have hit the onset of a segment
            
            section_hold = [section_hold;section{join}];
            spike_peaks_filt_hold = [spike_peaks_filt_hold spike_peaks_filt(join)];
            raw_peaks_hold = [raw_peaks_hold raw_peaks(join)];
            spike_peaks_hold = [spike_peaks_hold spike_peaks(join)];
            
        end

   if join+1<=length(find_segments_new)
        if find_segments_new(join)==0 | find_segments_new(join+1)==0 % If we aren't combining these sections, OR are coming off the end of combining sections
            
                new_events = new_events+1; % Update the new number of events
            
            if find_segments_new(join)==1 % If we are coming off a combined segment
                
                section_temp(new_events) = {section_hold};
                ratio_temp = ratio(cell2mat({section_hold}));
                temp_ind = find(ratio_temp==max(ratio_temp), 1); % MODIFIED: added 1
                temp_inds = cell2mat(section_temp(new_events));
                spike_peaks_temp(new_events) = temp_inds(temp_ind); % Gives the new ratio peak of the combined events
                
                % Now, search within window for the corresponding raw peak
                
                on_deck = data(spike_peaks_temp(new_events)-round(fs/search_width):spike_peaks_temp(new_events)+round(fs/search_width),tr);
                temp_ind = find(on_deck==max(on_deck));
                raw_peaks_temp(new_events) = temp_ind+spike_peaks_temp(new_events)-round(fs/search_width)-1;
                
                % Now, search within window for the corresponding filtered peak
                
                on_deck = smooth_wave(spike_peaks_temp(new_events)-round(fs/search_width):spike_peaks_temp(new_events)+round(fs/search_width),tr);
                temp_ind = find(on_deck==max(on_deck));
                spike_peaks_filt_temp(new_events) = temp_ind+spike_peaks_temp(new_events)-round(fs/search_width)-1;
                                
            else
                
                section_temp(new_events) = section(join);
                spike_peaks_filt_temp(new_events) = spike_peaks_filt(join);
                raw_peaks_temp(new_events) = raw_peaks(join);
                spike_peaks_temp(new_events) = spike_peaks(join);
                
            end
            
                section_hold = [];
                spike_peaks_filt_hold = [];
                raw_peaks_hold = [];
                spike_peaks_hold = [];               
                
        end
   end % end for ensuring we don't try and join the last segment to a segment that does not exist and throw an error
                           
    end % For join loop
    
    if new_events~=0 % If there were no sections to combine
                section = section_temp;
                spike_peaks_filt = spike_peaks_filt_temp;
                raw_peaks = raw_peaks_temp;
                spike_peaks = spike_peaks_temp;
    end
    
end % End for if we are applying min_dur
