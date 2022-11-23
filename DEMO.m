% Written by Stephen Evans swenceslaoe@hotmail.com
% % Copy pasted from the excel sheet '/Users/stephenevans/Desktop/Temp_dump_ASAP4_detection'

load_sample_data;

%% Chop Data and Normalize it with lowpass

chopped_ASAP3_1000hz = ASAP3_1000hz(1100:end);
chopped_ASAP3_500hz = ASAP3_500hz(550:end);

Fpass = 3;
fs = 1000;
[a,b] = butter(2,Fpass*(2*1/fs), 'low');
lowpass = filtfilt(a,b,chopped_ASAP3_1000hz); % filter the trace

Norm_and_chopped_ASAP3_1000hz = chopped_ASAP3_1000hz./lowpass;

Fpass = 3;
fs = 500;
[a,b] = butter(2,Fpass*(2*1/fs), 'low');
lowpass = filtfilt(a,b,chopped_ASAP3_500hz); % filter the trace

Norm_and_chopped_ASAP3_500hz = chopped_ASAP3_500hz./lowpass;
%figure;plot(chopped_ASAP3_500hz);hold on;plot(lowpass,'r');hold off

time = 0.002:0.002:length(chopped_ASAP3_500hz)/500;
points = 0.001:0.001:(length(time)-0.5)/500;

linear_int_500_to_1000hz = interp1(time,chopped_ASAP3_500hz,points,'linear');
linear_int_500_to_1000hz = linear_int_500_to_1000hz';

if sum(isnan(linear_int_500_to_1000hz(1)))==1
    linear_int_500_to_1000hz(1) = linear_int_500_to_1000hz(2);
end

Fpass = 3;
fs = 1000;
[a,b] = butter(2,Fpass*(2*1/fs), 'low');
lowpass = filtfilt(a,b,linear_int_500_to_1000hz); % filter the trace
%figure;plot(linear_int_500_to_1000hz);hold on;plot(lowpass,'r');hold off

Norm_linear_int_500_to_1000hz = linear_int_500_to_1000hz./lowpass;

t1 = 0.001:0.001:length(Norm_linear_int_500_to_1000hz)/1000;
t2 = 0.002:0.002:length(Norm_and_chopped_ASAP3_500hz)/500;
%figure;plot(t1,Norm_linear_int_500_to_1000hz,'k');hold on;plot(t2,Norm_and_chopped_ASAP3_500hz,'r');hold off;

%% Pull templates for 1000hz
% Here spike templates have been chosen from ASAP3, but the detection algorithm can take any arbitrary signal.

pre_spike = 50;
post_spike = 50;
sub_pre = 1;
sub_post = 12;
templates = [];
templates_hold = [];
templatess = [];
fs = 1000;
% [a,b] = butter(2,Fpass*(2*1/fs), 'high');

% From ASAP3_MM72 trial 1
template_peaks = [84,216,418,471,516,681,945,1001,1133,1201,1323,1405];

for i = 1:length(template_peaks)
    templates_hold(i,:) = Norm_and_chopped_ASAP3_1000hz(template_peaks(i)-pre_spike:template_peaks(i)+post_spike);
    %templatess(i,:) = filtfilt(a,b,templates_hold(i,:));% filter the trace
    templates(i,:) = templates_hold(i,pre_spike-sub_pre:post_spike+sub_post);
    %templates(i,:) = templates(i,:)+templates_hold(i,1);
    templates(i,:) = templates(i,:)./templates(i,1);
end

figure;plot(mean(templates))

%% Run detection with min_int = 10/fs

desired_FP = 1/(length(Norm_and_chopped_ASAP3_1000hz)*20);
cutoff = 0;
fs = 1000;
pre_spike = 50;
post_spike = 50;
Fpass = 20;

dbstop if error
[spikes_struct] = detect_events(templates,Norm_and_chopped_ASAP3_1000hz,desired_FP,pre_spike,post_spike,cutoff,fs,'Fpass',Fpass,'ratio_thresh',log(20*((length(Norm_and_chopped_ASAP3_1000hz)-(2*fs*cutoff))/size(templates,2))),'No_Filt','min_interval',10/fs);

% The threshold for the log likelihood ratio was chosen in a statistical
% manner, such that there is a 1/20 chance of getting a false positive
% every length(template) samples after accounting for the parts of the data
% vector that are not analyzed.
