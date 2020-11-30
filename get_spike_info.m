function[spikes_struct,whole_spike,Spike_Stack_norm,spikes,ratio_range,ratio_PDF,PDF] = get_spike_info(data,smooth_wave,ff,fflow,fF,ratio,fs,pre_spike,post_spike,spike_peaks,raw_peaks,spike_peaks_filt,W,section,updown,spikes_struct,tr,varargin);

nyquist_freq = fs/2;             % Nyquist frequency
cutoff_freqs = [4 11]; % Pull out theta band
Wn=cutoff_freqs/nyquist_freq;    % non-dimensional frequency
order = 1;
[filtb,filta]=butter(order,Wn,'bandpass'); % construct the filter

% Manually delete APs - Just ignore the part on the amplitude where spike
% is present. Didn't change analysis for findings any so dropped it.
% raw = [Mouse_2_120120_File_1(1:2170,2);Mouse_2_120120_File_1(2271:3911,2);Mouse_2_120120_File_1(4013:4318,2);Mouse_2_120120_File_1(4419:6479,2);Mouse_2_120120_File_1(6580:8532,2);Mouse_2_120120_File_1(8632:11121,2);Mouse_2_120120_File_1(11221:end,2)];

t_vec = 1/fs:1/fs:(pre_spike+post_spike+1)/fs;

raw_filt = []; % Filtered with butterworth bandpass at 4-11hz.
raw_hilbert = [];
phi = [];
phi_stack = [];
amp = [];
amp_stack = []; % Found by taking absolute value of hilbert transform
amp_stack_norm = [];
AP_height_norm = []; % Taken by finding max of raw data in section where spike is being taken.

    % bandpass filter, get hilbert transform
    raw = data(:,tr); % The raw vector
    raw_filt(tr,:)=filtfilt(filtb,filta,raw);   % filter the data with zero phase 
    % Hilbert
    raw_hilbert(tr,:) = hilbert(raw_filt(tr,:));
    t = 1/fs:1/fs:length(raw_hilbert(tr,:))/fs;
    % Compute phase of analytic signal
    phi(tr,:) = angle(raw_hilbert(tr,:));
    % Compute amplitude envelope
    amp(tr,:) = abs(raw_hilbert(tr,:));
                     
               % Now grab W around each spike               
               raw_spikes = [];
               peri = [];
               peri_raw = [];
               spikes = [];
               Spike_Stack_norm = []; 
               Spike_stack_norm_onset_aligned = [];
               Filtered_Spike_Stack_norm = []; 
               whole_spike = 0-pre_spike:1:post_spike;
               whole_spike_onset = 0-pre_spike:1:(post_spike*3);
               raw_norm = data(:,tr); %./fflow; % Get rid of slow drifts due to long plateus and photobleaching.
               amp_stack = [];
               amp_stack_norm = [];
               phi_stack = [];
               AP_height_norm = [];
               AP_height = [];
               
               spike_onset_times = [];
               onset_mat = [];
               onset_half_width = [];
               offset_half_width = [];
               offset_mat = [];
               total_width = [];
               total_half_width = [];
               raw_spike_integral = []; % Take sum of points over the "total" region, but in the raw trace
               filtered_spike_integral = []; % Take sum of points where data has been fitted for spike.
               
               % Get Theta Data new way (old way used same start and end inds).
                base_start = [];
                base_end = [];
                amp_base = [];
                amp_peak = []; 
                delta_theta = [];

                 for s = 1:length(spike_peaks_filt)                               % For each discovered spike
                    if spike_peaks(s)>W(end) && spike_peaks(s)<length(ff)+W(1) % Don't include weirdness at the start and end of the filtered trace                        
                        peri(s,:) = spike_peaks(s) + W;                          % Grab points around spike peak
                        spikes(s,:) = fF(peri(s,:),tr);                        % Store the perispike filtered trace                        
                        peri_raw(s,:) = raw_peaks(s) + W;
                        raw_spikes(s,:) = data(peri_raw(s,:),tr);
                        
                        waveform(s,:) = raw_peaks(s) + whole_spike;
                        ratio_waveform(s,:) = spike_peaks(s) + whole_spike;
                        
                        if waveform(s,end)>length(raw_norm) % If the waveform goes over the edge of the vector
                            new_temp = raw_norm(waveform(s,1):end)';new_temp2 = ones((waveform(s,end)-length(raw_norm)),1)';                            
                            Spike_Stack_norm(s,:) = [new_temp new_temp2]; 
                            new_temp = data(waveform(s,1):end,tr)';new_temp2 = ones((waveform(s,end)-length(raw_norm)),1)';
                            Spike_Stack_original(s,:) = [new_temp new_temp2];                             
                            new_temp = ratio(waveform(s,1):end)';new_temp2 = ones((waveform(s,end)-length(ratio)),1)'; 
                            ratio_spikes(s,:) = [new_temp new_temp2];                             
                            
                        else Spike_Stack_norm(s,:) = raw_norm(waveform(s,:))./mean(raw_norm(waveform(s,1))); 
                            Spike_Stack_original(s,:) = data(waveform(s,:),tr);  
                            ratio_spikes(s,:) = ratio(ratio_waveform(s,:),tr);
                        end                        
                        
                            waveform_filt(s,:) = spike_peaks_filt(s) + whole_spike;
                            
                            if waveform_filt(s,end)>length(amp) % If the waveform goes over the edge of the vector
                                new_temp = smooth_wave(waveform_filt(s,1):end,tr)'; new_temp2 = ones((waveform_filt(s,end)-length(smooth_wave)),1)';
                                Filtered_Spike_Stack_norm(s,:) = [new_temp new_temp2];
                                new_temp = amp(tr,waveform_filt(s,1):end); new_temp2 = ones((waveform_filt(s,end)-length(smooth_wave)),1)';
                                amp_stack(s,:) = [new_temp new_temp2];
                                amp_stack_norm(s,:) = amp_stack(s,:)./mean(amp_stack(s,:));
                                new_temp = phi(tr,waveform_filt(s,1):end); new_temp2 = ones((waveform_filt(s,end)-length(smooth_wave)),1)';
                                phi_stack(s,:) = [new_temp new_temp2];
                            else
                                Filtered_Spike_Stack_norm(s,:) = smooth_wave(waveform_filt(s,:),tr); 
                                amp_stack(s,:) = amp(tr,waveform(s,:));
                                amp_stack_norm(s,:) = amp_stack(s,:)./mean(amp_stack(s,:));
                                phi_stack(s,:) = phi(tr,waveform(s,:));
                            end
                            
           % Begin width section
           % Find height of AP of nonfiltered spike
           if sum(strcmp(varargin{1},'No_Filt'))>0 % If the data and templates are optimally filtered
                AP_height(s,:) = data(raw_peaks(s),tr); 
                peak_ind = pre_spike+1;
                half_height = (AP_height(s,:)+1)/2;
           else peak_ind = pre_spike+1; AP_height_norm(s,:) = Spike_Stack_norm(s,peak_ind);half_height = (AP_height_norm(s,:)+1)/2;
           end
                                                                                                                                                
                        % calculate by Ramp amplitude was quantified as the difference between peak 
                        % (10 cm average around most depolarized value) and the baseline 
                        % (10 cm average around most hyperpolarized value).
                    
                        % The below while loops move backwards and forwards from the peak until
                        % the trace stops going down. This is called the total_width. Halfway 
                        % between the peak and 1 is called the half_width.   
                        
                        % Move back and forth in ratio
                        
                        onset_advance = 1; % enter while loop                       
                        ratio_onset_count = 1;
                        
                        while onset_advance==1 % Find total spike width by figuring out how far until it flips direction in filtered spikes

                                onset_advance = ratio_spikes(s,peak_ind-ratio_onset_count)-ratio_spikes(s,peak_ind-(ratio_onset_count-1))<0;

                            ratio_onset_count = ratio_onset_count+1;
                            
                            if onset_advance==0
                                ratio_onset_count = ratio_onset_count-2; % Correct for the counter and the way it differences
                            end
                            
                            if peak_ind-ratio_onset_count==0; onset_advance = 0;ratio_onset_count = ratio_onset_count-1; end
                        end
                        
                        ratio_onset_mat(s) = ratio_onset_count*(1/fs); % keep track of ms to peak
                        
                        offset_advance = 1; % enter while loop
                        ratio_offset_count = 1;
                        while offset_advance==1
                            
                            if peak_ind+ratio_offset_count<=size(ratio_spikes(s,:),2) % As long as we aren't out of bounds
                                offset_advance = ratio_spikes(s,(peak_ind+ratio_offset_count))-ratio_spikes(s,peak_ind+(ratio_offset_count-1))<0;
                                ratio_offset_count = ratio_offset_count+1;
                                if offset_advance==0
                                    ratio_offset_count = ratio_offset_count-2; % Correct for counter and the way it does differencing
                                end
                            end
                            if peak_ind+ratio_offset_count>size(ratio_spikes(s,:),2);offset_advance = 0;
                                ratio_offset_count = ratio_offset_count-1;
                            end
                        end
                        
                        ratio_offset_mat(s) = ratio_offset_count*(1/fs); % keep track of ms to peak
                        % Move back and forth on filtered data

                        onset_advance = 1; % enter while loop                       
                        onset_count = 1;
                        
                        while onset_advance==1 % Find total spike width by figuring out how far until it flips direction in filtered spikes
                            if updown==1 % If the event goes up
                                onset_advance = Filtered_Spike_Stack_norm(s,peak_ind-onset_count)-Filtered_Spike_Stack_norm(s,peak_ind-(onset_count-1))<0;
                            end
                            if updown==0
                                onset_advance = Filtered_Spike_Stack_norm(s,peak_ind-onset_count)-Filtered_Spike_Stack_norm(s,peak_ind-(onset_count-1))>0;
                            end
                            onset_count = onset_count+1;
                            
                            if onset_advance==0
                                onset_count = onset_count-2; % Correct for the counter and the way it differences
                            end
                            
                            if peak_ind-onset_count==0; onset_advance = 0;onset_count = onset_count-1; end
                        end
                        
                        onset_mat(s) = onset_count*(1/fs); % keep track of ms to peak
  
                        offset_advance = 1; % enter while loop
                        offset_count = 1;
                        while offset_advance==1
                            
                            if peak_ind+offset_count<=size(Filtered_Spike_Stack_norm(s,:),2) % As long as we aren't out of bounds
                                offset_advance = Filtered_Spike_Stack_norm(s,(peak_ind+offset_count))-Filtered_Spike_Stack_norm(s,peak_ind+(offset_count-1))<0;
                                offset_count = offset_count+1;
                                if offset_advance==0
                                    offset_count = offset_count-2; % Correct for counter and the way it does differencing
                                end
                            end
                            if peak_ind+offset_count>size(Filtered_Spike_Stack_norm(s,:),2);offset_advance = 0;
                                offset_count = offset_count-1;
                            end
                        end

                        [~,half_ind] = min(abs(Spike_Stack_original(s,(peak_ind-round(fs/50)):peak_ind)-half_height)); % Find the point closest to half height.
                        half_ind_on = half_ind + (peak_ind - round(fs/50));
                        onset_half_width(s) = (peak_ind - half_ind_on)*(1/fs); % number of points back from peak it hits closest to half. Resolution limited by sampling rate.  

                        [~,half_ind] = min(abs(Spike_Stack_original(s,(peak_ind:peak_ind+round(fs/50)))-half_height)); % Find the point closest to half height.
                        half_ind_off = half_ind + peak_ind;
                        offset_half_width(s) = (half_ind_off - peak_ind)*(1/fs); % number of points back from peak it hits closest to half. Resolution limited by sampling rate.                                                                                                                                                                 

                        offset_mat(s) = offset_count*(1/fs); % keep track of ms to peak
                        total_width(s) = (onset_mat(s) + offset_mat(s))*1000; % In ms
                        total_half_width(s) = (onset_half_width(s)+offset_half_width(s))*1000; % In ms 
                        ratio_full_width(s) = ratio_onset_mat(s) + ratio_offset_mat(s) + 1/fs;;
                        ratio_onset(s) = spike_peaks(s) - (ratio_onset_mat(s)*fs);
                        filtered_spike_integral(s) = sum(smooth_wave(spike_peaks_filt(s)-(onset_mat(s)*fs):spike_peaks_filt(s)+(offset_mat(s)*fs),tr)) - ((onset_mat(s)*fs)+(offset_mat(s)*fs)+1);
                        % Subtracting the length of the vector assumes it
                        % is normalized to have background value of 1!
                        raw_spike_integral(s) = sum(data(raw_peaks(s)-(onset_mat(s)*fs):raw_peaks(s)+(offset_mat(s)*fs),tr)) - ((onset_mat(s)*fs)+(offset_mat(s)*fs)+1);
                        % Subtracting the length of the vector assumes it
                        % is normalized to have background value of 1!
                        % End width section
                        
                        [~,I] = min(amp_stack_norm(s,round(0.25*pre_spike):pre_spike)); % Find theta min before spike.
                        base_start(s) = I+round(0.25*pre_spike)-1; % Shift ind to trace value
                        base_end(s) = base_start(s)+(round(fs*(2/317)));
                        time_before_spike_theta(s,1) = ((pre_spike-(base_start(s)-1))/fs)*1000; % In ms

                        amp_base(s) = mean(amp_stack_norm(s,base_start(s):base_end(s)));
                        amp_peak(s) = mean(amp_stack_norm(s,pre_spike-round(fs*(2/317)))); % Picks a spot 9ms before peak.
                        delta_theta(s) = amp_peak(s)-amp_base(s);
                        
                    end % end for if statement that excludes trace ends
                    
                end   % end for s loop (spikes)            
                              
                
                ratio_range = min(ratio(ratio(:,tr)~=-Inf,tr)):max(ratio(ratio(:,tr)~=Inf,tr));                     
                ratio_PDF = ratio(ratio(:,tr)~=-Inf,tr);
                ratio_PDF = ratio_PDF(ratio_PDF~=Inf);
                PDF = histcounts(ratio_PDF, ratio_range);                       % Get probclcability density function
                PDF = PDF/sum(PDF); 

   if length(spike_peaks)>0
       
        delta_theta_norm = delta_theta-mean(delta_theta);
        delta_theta_norm = delta_theta_norm./std(delta_theta_norm);

        % Get Voltage (Fluorescence) Ramp Data

        ramp_start = pre_spike-round(fs*(64/317));if ramp_start<0 ramp_start = 1;end
        ramp_end = ramp_start+round(fs*(3/317));
        ramp_base = mean(Spike_Stack_norm(:,ramp_start:ramp_end),2);
        ramp_peak = mean(Spike_Stack_norm(:,pre_spike-round(fs*(6/317)):pre_spike-round(fs*(3/317))),2);
        delta_ramp = ramp_peak-ramp_base; % This is the change between a 25 ms average from 205 to 180 ms
        % before the spike to the average voltage 35 to 10 ms before the spike.

        delta_ramp_norm = delta_ramp-mean(delta_ramp);
        delta_ramp_norm = delta_ramp_norm./std(delta_ramp_norm);

        ramp_theta_interaction = delta_theta_norm'.*delta_ramp_norm;

        theta_ramp_mat = [delta_theta_norm',delta_ramp_norm,ramp_theta_interaction,phi_stack(:,pre_spike)];
        % Note that p values change on the terms depending on what terms are
        % included!
   end
   
           if isempty(spike_peaks)==0 % If we found spikes 
            spikes_struct(tr).ratio_full_width = ratio_full_width; % Full width of the spike in the ratio vector  
            spikes_struct(tr).ratio_based_onset_detection = ratio_onset;
            spikes_struct(tr).filtered_spike_integral = filtered_spike_integral; 
            spikes_struct(tr).raw_spike_integral = raw_spike_integral;     
            spikes_struct(tr).AP_height = AP_height;   
            spikes_struct(tr).delta_theta = delta_theta;
            spikes_struct(tr).spikes = spikes; % The spike traces pulled out from the highpass filtered data
            spikes_struct(tr).raw_spikes = raw_spikes; % The spike traces pulled out from the raw data
            spikes_struct(tr).raw_spike_peaks = raw_peaks; % From the filtered raw trace
            spikes_struct(tr).spike_peaks = spike_peaks; % From the peaks in ratio
            spikes_struct(tr).raw_spike_inds = peri_raw;
            spikes_struct(tr).spike_inds = peri; % The inds for the spikes. From the highpass filtered
            % data without the ends chopped, so corresponds to the real data.
            spikes_struct(tr).above_threshold = section; % This contains the inds in the raw data
            % for the data points that are above ratio_thresh, without any end
            % chopping. Note that there are usually a few fake spikes in the
            % first 0.1 seconds from the filtering artifact, so you should
            % manually drop any values that are less than 0.1*fs.
            
            spikes_struct(tr).Full_Spike_Traces_norm = Spike_Stack_norm;
            spikes_struct(tr).AP_height_norm = AP_height_norm;
            spikes_struct(tr).Filtered_Full_Spike_Traces_norm = Filtered_Spike_Stack_norm;
            spikes_struct(tr).total_width = total_width;
            spikes_struct(tr).total_half_width = total_half_width;
            spikes_struct(tr).delta_theta_norm = delta_theta_norm;
            spikes_struct(tr).fluorescence_ramp_start = ramp_start;
            spikes_struct(tr).fluorescence_ramp_end = ramp_end;
            spikes_struct(tr).delta_fluorescence_ramp = delta_ramp;
            spikes_struct(tr).fluorescence_ramp_base = ramp_base;
            spikes_struct(tr).fluorescence_ramp_peak = ramp_peak;
            spikes_struct(tr).delta_ramp_norm = delta_ramp_norm;
            spikes_struct(tr).ramp_theta_interaction = ramp_theta_interaction;
            spikes_struct(tr).theta_ramp_mat = theta_ramp_mat;
            spikes_struct(tr).amp_stack = amp_stack;
            spikes_struct(tr).amp_stack_norm = amp_stack_norm;
            spikes_struct(tr).phi_stack = phi_stack;
            spikes_struct(tr).onset_mat = onset_mat;
            spikes_struct(tr).onset_half_width = onset_half_width;
            spikes_struct(tr).offset_half_width = offset_half_width;
            spikes_struct(tr).offset_mat = offset_mat;                      
        end