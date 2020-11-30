%% This section looks at the data stored in Google drive --> Shared drives --> Stephen Evans shared drive --> ASAP4_paper --> Adrians_data --> cell AN080819-1
%% AN080819-1 was found to have nothing, and according to Adrian is trash. Section deleted.
%% AN070120 is ASAP4.4, and has lots of events
%% For AN070120, ASAP4.4 data
% Volumes --> Mac Backup --> Adrians_data -->
% AN070120_selection-export-signals_ASAP4.4

% Load the data

filelist = dir('*.mat');
files = struct;

for i = 1:length(filelist)
    
   files(i).file = load([filelist(i).name]);
    
end

%% Look at the data

%1)	2?nd order butterworth highpass at 0.5hz. 
%2)	Downsample from 440 to 220 hz by summing.
%3)	11 point moving average

data = [];

for i = 1:length(files)
   
    fs = files(i).file.imaging_fs;
    [a,b] = butter(2,0.5*(2*1/fs), 'high'); % get rid of slow drift
            % GC                     
            gc_raw = files(i).file.imaging_ROIs_dendrite_d1_Ch1_sig;% - gc_back_fit;
            gc_raw = gc_raw(logical(abs((isnan(gc_raw)-1))));
            t = 1:1:length(gc_raw);  

            [s] = exp2fit(t,double(gc_raw),1);
            gc_fit = s(1) + s(2)*exp(-t/s(3));
            gc = gc_raw./gc_fit; % To get in terms of df/F
            gc = filtfilt(a,b,double(gc)); % To get rid of slow drift
            gc = gc+1; % To get back to df/F
            
            first = downsample(gc,2);
            second = downsample(gc(2:end),2);
        
            if length(first)>length(second)
                first = first(1:length(second));
            end
            if length(second)>length(first)
                second = second(1:length(first));
            end

            gc = first+second; % To clean up some noise
            gc = gc./mean(gc); % To get it back to df/F
            
            % RC           
            rc_raw = files(i).file.imaging_ROIs_dendrite_d1_Ch2_sig;% - rc_back_fit;
            rc_raw = rc_raw(logical(abs((isnan(rc_raw)-1))));
            t = 1:1:length(rc_raw); 
            
            [s] = exp2fit(t,double(rc_raw),1);
            rc_fit = s(1) + s(2)*exp(-t/s(3));
            rc = rc_raw./rc_fit; % To get in terms of df/F
            rc = filtfilt(a,b,double(rc)); % To get rid of slow drift
            rc = rc+1; % To get back to df/F
            
            first = downsample(rc,2);
            second = downsample(rc(2:end),2);
        
            if length(first)>length(second)
                first = first(1:length(second));
            end
            if length(second)>length(first)
                second = second(1:length(first));
            end

            rc = first+second; % To clean up some noise
            rc = rc./mean(rc); % To get it back to df/F
            ratio = gc./rc;
            ratio = evans_conv([1,1,1,1,1,1,1,1,1,1,1],ratio,1,0);
            
        %figure;plot(ratio)
            time = 1/(files(1).file.imaging_fs/2):1/(files(1).file.imaging_fs/2):(length(gc)/1/(files(1).file.imaging_fs/2));
            figure;title('ASAP4.4');hold on;         
            subplot(3,1,1);plot(time(100:end-200),gc(100:end-200),'g');title('Green Channel');         
            subplot(3,1,2);plot(time(100:end-200),rc(100:end-200),'r');title('Red Channel');       
            subplot(3,1,3);plot(time(100:end-200),ratio(100:end-200),'k');title('Ratio');xlabel('Time Seconds');ylabel('Filtered and Normalized')  

        data(:,i) = ratio(1:9894)';   
    
end

ASAP44_files = files;

%% Pick templates out from ratio, run detection

pre_spike = 20;
post_spike = 15;

templates1 = [];
templates2 = [];
templates = [];

template_peaks1 = [1247,5679,6248,9228,9639];

for i = 1:length(template_peaks1)
    templates1(i,:) = data(template_peaks1(i)-pre_spike:template_peaks1(i)+post_spike,1);
    templates1(i,:) = templates1(i,:)./templates1(i,1);
end

template_peaks2 = [1201,1960,4055,5256,5685,9752];
for i = 1:length(template_peaks2)
    templates2(i,:) = data(template_peaks2(i)-pre_spike:template_peaks2(i)+post_spike,2);
    templates2(i,:) = templates2(i,:)./templates2(i,1);
end

templates = [templates1;templates2];
figure;plot(mean(templates)) % looks great!

%% Run Detection - Amazing!

dbstop if error
desired_FP = 1/(size(data,1)*20);
cutoff = 0.1;
pre_spike = 150;
post_spike = 90;
fs = files(i).file.imaging_fs/2; % From downsampling
thresh = log(20*(size(data,1)/size(templates,2)));
[spikes_struct] = detect_events(templates,data,desired_FP,pre_spike,post_spike,cutoff,fs,'Fpass',10,'No_Filt','ratio_thresh',thresh);

%% Look at spike integrals (need to clean out artifacts too). Looks perfectly poisson. Damn.

filtered_integral = [];
raw_integral = [];
for i = 1:length(spikes_struct)
   
    filtered_integral = [filtered_integral spikes_struct(i).filtered_spike_integral];
    raw_integral = [raw_integral spikes_struct(i).raw_spike_integral];
    
end
figure;hist(filtered_integral,[0.1:0.1:10])

%% For AN091619-1, ASAP3 data
% Volumes --> Mac Backup --> Adrians_data -->
% AN091619-1 selection-export-signals_ASAP3

% Load the data

filelist = dir('*.mat');
files = struct;

for i = 1:length(filelist)
    
   files(i).file = load([filelist(i).name]);
    
end

%% Look at the data

%1)	2'nd order butterworth highpass at 0.5hz. 
%2)	Downsample from 440 to 220 hz by summing.
%3)	11 point moving average

data = [];

for i = 1:length(files)
   
    fs = files(i).file.imaging_fs;
    [a,b] = butter(2,0.5*(2*1/fs), 'high'); % get rid of slow drift
            % GC                     
            gc_raw = files(i).file.imaging_ROIs_dendrite_d1_Ch1_sig;% - gc_back_fit;
            gc_raw = gc_raw(logical(abs((isnan(gc_raw)-1))));
            t = 1:1:length(gc_raw);  

            [s] = exp2fit(t,double(gc_raw),1);
            gc_fit = s(1) + s(2)*exp(-t/s(3));
            gc = gc_raw./gc_fit; % To get in terms of df/F
            gc = filtfilt(a,b,double(gc)); % To get rid of slow drift
            gc = gc+1; % To get back to df/F
            
            
            % RC           
            rc_raw = files(i).file.imaging_ROIs_dendrite_d1_Ch2_sig;% - rc_back_fit;
            rc_raw = rc_raw(logical(abs((isnan(rc_raw)-1))));
            t = 1:1:length(rc_raw); 
            
            [s] = exp2fit(t,double(rc_raw),1);
            rc_fit = s(1) + s(2)*exp(-t/s(3));
            rc = rc_raw./rc_fit; % To get in terms of df/F
            rc = filtfilt(a,b,double(rc)); % To get rid of slow drift
            rc = rc+1; % To get back to df/F
            
            
            ratio = gc./rc;
            ratio = evans_conv([1,1,1,1,1],ratio,1,0);
            
        %figure;plot(ratio)
        
            time = 1/(files(1).file.imaging_fs):1/(files(1).file.imaging_fs):(length(gc)/1/(files(1).file.imaging_fs));
            figure;title('ASAP3');hold on;         
            subplot(3,1,1);plot(time(100:end-200),gc(100:end-200),'g');title('Green Channel');         
            subplot(3,1,2);plot(time(100:end-200),rc(100:end-200),'r');title('Red Channel');       
            subplot(3,1,3);plot(time(100:end-200),ratio(100:end-200),'k');title('Ratio');xlabel('Time Seconds');ylabel('Filtered and Normalized')  


        data(:,i) = ratio(1:6990)';   
    
end

ASAP3_files = files;

%% Pick templates out from ratio, run detection

pre_spike = 13;
post_spike = 13;

templates1 = [];
templates2 = [];
templates = [];

template_peaks1 = [238,1084,1856];

for i = 1:length(template_peaks1)
    templates1(i,:) = data(template_peaks1(i)-pre_spike:template_peaks1(i)+post_spike,1);
    templates1(i,:) = templates1(i,:)./templates1(i,1);
end

template_peaks2 = [556,2772,3606,4314,4958];
for i = 1:length(template_peaks2)
    templates2(i,:) = data(template_peaks2(i)-pre_spike:template_peaks2(i)+post_spike,2);
    templates2(i,:) = templates2(i,:)./templates2(i,1);
end

template_peaks3 = [1182,1840,3277];
for i = 1:length(template_peaks3)
    templates3(i,:) = data(template_peaks3(i)-pre_spike:template_peaks3(i)+post_spike,3);
    templates3(i,:) = templates3(i,:)./templates3(i,1);
end

templates = [templates1;templates2;templates3];
figure;plot(mean(templates)) % looks great!

%% Run Detection - Amazing!

dbstop if error
desired_FP = 1/(size(data,1)*20);
cutoff = 0.1;
pre_spike = 150;
post_spike = 90;
thresh = log(20*(size(data,1)/size(templates,2)));
[spikes_struct] = detect_events(templates,data,desired_FP,pre_spike,post_spike,cutoff,fs,'Fpass',10,'No_Filt','ratio_thresh',thresh);

%% Pull out average waveform with std shading, half width and EPSP height histograms
% Don't forget to ignore events that are too close together, take one with
% biggest AP height.

% Load spikes_struct from Mac_Backup --> Adrians_data -->
% AN070120_selection-export-signals_ASAP4.4

ASAP44_half_widths = [];
ASAP44_total_widths = [];
ASAP44_spike_trace = [];
ASAP44_AP_height = [];

for i = 1:length(ASAP44_spikes_struct)
   
    ASAP44_half_widths = [ASAP44_half_widths ASAP44_spikes_struct(i).total_half_width];
    ASAP44_total_widths = [ASAP44_total_widths ASAP44_spikes_struct(i).total_width];
    ASAP44_AP_height = [ASAP44_AP_height;ASAP44_spikes_struct(i).AP_height];
    ASAP44_spike_trace = [ASAP44_spike_trace;ASAP44_spikes_struct(i).Full_Spike_Traces_norm];
    
end
fs = 220;
time = 1/fs:1/fs:length(ASAP44_spike_trace)/fs;
figure;plot(time,mean(ASAP44_spike_trace))
%stdshade(ASAP44_spike_trace)

%% Don't forget to ignore events that are too close together, take one with
% biggest AP height.

ASAP3_half_widths = [];
ASAP3_total_widths = [];
ASAP3_AP_height = [];
ASAP3_spike_trace = [];
ASAP3_height_percent = [];

for i = 1:length(ASAP3_spikes_struct)
   
    ASAP3_half_widths = [ASAP3_half_widths ASAP3_spikes_struct(i).total_half_width];
    ASAP3_total_widths = [ASAP3_total_widths ASAP3_spikes_struct(i).total_width];
    ASAP3_AP_height = [ASAP3_AP_height;ASAP3_spikes_struct(i).AP_height];
    ASAP3_height_percent = [ASAP3_height_percent;1-ASAP3_AP_height];
    ASAP3_spike_trace = [ASAP3_spike_trace;ASAP3_spikes_struct(i).Full_Spike_Traces_norm];
    
end
fs = 233;
time = 1/fs:1/fs:length(ASAP3_spike_trace)/fs;
figure;plot(time,mean(ASAP3_spike_trace))
%stdshade(ASAP3_spike_trace)

%% Make plots for meeting - total width histograms

figure;
Hist_ASAP44_half = histogram(ASAP44_total_widths);hold on;
Hist_ASAP3_half = histogram(ASAP3_total_widths);
box off
legend('ASAP44','ASAP3','location','northeast')
legend boxoff
title('Total Event Widths Histogram')

%% Half width histograms

figure;
Hist_ASAP44_half = histogram(ASAP44_half_widths);hold on;
Hist_ASAP3_half = histogram(ASAP3_half_widths);
box off
legend('ASAP44','ASAP3','location','northeast')
legend boxoff
title('Event Half Widths Histogram')

%% AP height histograms

figure;
Hist_ASAP44_AP_height = histogram(ASAP44_AP_height-1);hold on;
Hist_ASAP3_AP_height = histogram(1-ASAP3_AP_height);
box off
legend('ASAP44','ASAP3','location','northeast')
legend boxoff
title('Event Percent Fluorescence Change Histogram')

%% Show example traces ASAP3

time = 1/ASAP3_files(1).file.imaging_fs:1/ASAP3_files(1).file.imaging_fs:length(ASAP3_spike_trace)/ASAP3_files(1).file.imaging_fs;
stdshade(ASAP3_spike_trace,0.5,'r',time)
title('ASAP3 Trace and std')
ylabel('Normalized Fluorescence')
xlabel('Time seconds')

%% Show example traces ASAP44

time = 1/ASAP44_files(2).file.imaging_fs:1/ASAP44_files(2).file.imaging_fs:length(ASAP44_spike_trace)/ASAP44_files(2).file.imaging_fs;
stdshade(ASAP44_spike_trace,0.5,'b',time)
title('ASAP44 Trace and std')
ylabel('Normalized Fluorescence')
xlabel('Time seconds')

%% Show flipped ASAP3 trace with ASAP44 trace

ASAP3_time = 1/ASAP3_files(1).file.imaging_fs:1/ASAP3_files(1).file.imaging_fs:length(ASAP3_spike_trace)/ASAP3_files(1).file.imaging_fs;
stdshade((1-ASAP3_spike_trace),0.5,'r',ASAP3_time);hold on;
ASAP44_time = 1/(ASAP44_files(2).file.imaging_fs/2):1/(ASAP44_files(2).file.imaging_fs/2):length(ASAP44_spike_trace)/(ASAP44_files(2).file.imaging_fs/2);
ASAP44_time = ASAP44_time;
offset = ASAP44_time(151) - ASAP3_time(151);
stdshade(ASAP44_spike_trace-1,0.5,'b',ASAP44_time-offset)
title('ASAP3 Trace and std')
ylabel('Normalized Fluorescence')
xlabel('Time seconds')





