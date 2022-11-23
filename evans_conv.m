function[output] = evans_conv(filter,signal,dim,col_too)
% In the form [output] = evans_conv(filter,signal,dim,col_too);
% This function takes each element of filter, and multiplies corresponding
% elements in signal by them, then takes the mean of the result, creating a new point in
% output, before sliding over 1 element and repeating the process.
% The beginning and end of signal is zero padded. This generally assumes
% that filter is significantly smaller than signal. For example, suppose
% you have a 2 element filter, and a 1000 element signal. In the first
% iteration, the first element of filter will be multiples by zero, and the
% second by the first element of signal, then summed and divided by filter length, to give the first
% element of output. The filter then slides over 1 element, and now the
% first element of filter is multiplies by the first element of signal, the
% second element of filter by the second of signal, the result summed and
% divided by 2 and
% turned into the second element of output. The process repeats until the
% end of signal is reached. We're just taking the dot product. The answer
% is then divided by the length of the filter.

% If dim is 2, the function slides the filter across each row, then each
% column, then adds them, if the variable col_too is a 1. Otherwise it just
% does rows.

% Currently sums matrices, then turns into a logical matrix.

if dim==1

if size(filter,1)~=1  % Make sure dimensions are the same
    filter = filter';
end

if size(signal,1)~=1 % Make sure dimensions are the same
    signal = signal';
end

zero_pad = zeros(1,length(filter)-1); % Number of zero elements to pad signal by.
new_signal = [zero_pad signal zero_pad];

output = []; % Vector to be filled

for i = 1:length(signal)

    transposed = new_signal(i:i+length(filter)-1)';
    output(i) =  filter*transposed; % Take the dot product.
end
end

if dim ==2
        if size(filter,1)~=1  % Make sure dimensions are the same
            filter = filter';
        end
    for i = 1:size(signal,1) % Count rows
        
        on_deck = signal(i,:); % Take i'th row of signal matrix
        zero_pad = zeros(1,length(filter)-1); % Number of zero elements to pad signal by.
        new_signal = [zero_pad on_deck zero_pad];
        
        for ii = 1:length(on_deck) % Convolve
            transposed = new_signal(ii:ii+length(filter)-1)';
            output_rows(i,ii) =  filter*double(transposed); % Take the dot product.
        end
        
    end % End for i loop
    
if col_too==1
% Now for columns

    for i = 1:size(signal,2) % Count columns
        
        on_deck = signal(:,i); % Take i'th row of signal matrix
        on_deck = on_deck';
        zero_pad = zeros(1,length(filter)-1); % Number of zero elements to pad signal by.
        new_signal = [zero_pad on_deck zero_pad];
        
        for ii = 1:length(on_deck) % Convolve
            transposed = new_signal(ii:ii+length(filter)-1)';
            output_columns(i,ii) =  filter*transposed; % Take the dot product.
        end
        
    end % End for i loop
    
    output = output_columns'+output_rows; % Sum the matrices
    output = abs(output); output = output>0;% Turn into logical matrix
end % end for if col_too==1
    
    if col_too~=1
        output = output_rows;
    end

end % End for if dim==2

output = output./length(filter);
