function baseline = baseline4peaks(x, fs)

% iterative lowpass filtering for traces with many peaks

F_baseline = 0.3; % baseline pass freq 0.3
[a,b] = butter(2,F_baseline*(2*1/fs), 'low');
Flp = filtfilt(a,b,x); % filter the trace
F = x./Flp;

mask = F > max(threshold(F), 1+2*EMstd(F));

F = x;
F(mask) = mean(x(~mask));

F_baseline = 0.5; % baseline pass freq 0.5
[a,b] = butter(2,F_baseline*(2*1/fs), 'low');
baseline = filtfilt(a,b,F); % filter the trace