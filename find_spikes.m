function[spike_peaks_filt,raw_peaks,spike_peaks,section,closest_to_start] = find_spikes(ASAPtime,cutoff,ratio,ff,fF,smooth_wave,data,ratio_thresh,updown,pre_spike,post_spike,fs,tr);

[~,closest_to_start] = min(abs(ASAPtime-cutoff));
                temp = find(ratio(:,tr)>ratio_thresh);
                Above_thresh = temp(temp>=closest_to_start);
                section = {};
                spike_peaks = [];
                spikes = [];
                raw_spikes = [];
                ratio_range = [];
                spike_peaks_filt = [];
                raw_peaks = [];

                if isempty(Above_thresh)==0 % If there are spikes detected
                    
                    diff_above_thresh = diff(Above_thresh);
                    logical_above_thresh = diff_above_thresh==1; % Clusters
                    % of ones correspond to a single spiking event, noted
                    % by the corresponding inds in Above_thresh. Note that
                    % the ind after the last 1 of each cluster is also
                    % included in the previous cluster, because of the way
                    % differencing works.

                    ratio_spike_peak_inds = Above_thresh(diff_above_thresh>1); % Using the point where ratio first crosses threshold to count as spike onset.
                    % Not great because it will vary dramatically depending
                    % on how the threshold is set.
                    not_zero = find(logical_above_thresh~=0); % Gives indices for Nonzero Elements
                    if isempty(not_zero)==0
                        start_inds_above_thresh = unique([not_zero(1);not_zero(diff([0;not_zero])>1)]); % Segment Start Indices
                    else start_inds_above_thresh = [];
                    end
                        Above_thresh_duplicate = Above_thresh; % Will use this to catch cases where there are multiple
                    % zeros in a row in logical_above_thresh, which
                    % represents multiple points where the log likelihood
                    % only crossed threshold for a single data point.
                    
                    for k1 = 1:length(start_inds_above_thresh) % 1 loop for each start detected
                       
                        on_deck = start_inds_above_thresh(k1);
                        counter = 0;

                            while logical_above_thresh(on_deck+counter)==1
                                counter = counter+1;
                                if on_deck+counter>=length(logical_above_thresh);
                                    break;
                                end % get out of the while loop. logical
                                % above_thresh ends on a 1.
                            end
                            end_inds_above_thresh(k1) = on_deck+counter; % Segment End Indices 
                            section{k1} = [Above_thresh(start_inds_above_thresh(k1):end_inds_above_thresh(k1))];
                            Above_thresh_duplicate(start_inds_above_thresh(k1):end_inds_above_thresh(k1)) = 0;
                            [~,spike_peaks(k1)] = max(ratio(section{k1},tr)); % Find ind of max within spiking event
                            spike_peaks(k1) = spike_peaks(k1)+Above_thresh(start_inds_above_thresh(k1))-1; % Adjust ind
                            % to fit into original vector..
                    end
                    if isempty(k1)==1
                        k1 = 0;
                    end
                spike_peaks(k1+1:k1+sum(Above_thresh_duplicate>0)) = Above_thresh_duplicate(Above_thresh_duplicate>0);
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
                      search_width = 15; % Smaller = search a wider area for peak
                      
                     for j = 1:length(spike_peaks) % 1 loop for each spike. Peak in trace isn't always the same
                         % as the peak in ratio, sometimes off by a few
                         % sample points, so we have to shift it.
                         
                        on_deck = data(spike_peaks(j)-round(fs/search_width):spike_peaks(j)+round(fs/search_width),tr);
                        
                        if updown == 1 % Depends on if the indicator goes up or down
                            [raw_pks(j), raw_locs(j)] = max(on_deck);
                        else
                            [raw_pks(j), raw_locs(j)] = min(on_deck);
                        end
                        raw_peaks_temp = spike_peaks(j)-(round(fs/search_width)+1)+raw_locs(j);
                        if raw_peaks_temp-pre_spike>0
                            spike_count_raw = spike_count_raw+1;
                            raw_peaks(spike_count_raw) = raw_peaks_temp;                       
                        end
                        
                        on_deck_filt = smooth_wave(spike_peaks(j)-round(fs/search_width):spike_peaks(j)+round(fs/search_width),tr);
                        if updown==1
                            [pks(j), locs(j)] = max(on_deck_filt);
                        else
                            [pks(j), locs(j)] = min(on_deck_filt);
                        end
                        
                        spike_peaks_filt_temp = spike_peaks(j)-(round(fs/search_width)+1)+locs(j);
                        if spike_peaks_filt_temp-pre_spike>0
                            spike_count_filt = spike_count_filt+1;
                            spike_peaks_filt(spike_count_filt) = spike_peaks_filt_temp;
                        end
                        
                     end % end for j loop
                end % end for if we found spikes