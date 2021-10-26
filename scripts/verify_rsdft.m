% Verify the functionality of the implementation of the rSDFT
% This is done by comparing the generated output from the HW model with
% the output from the MATLAB model.

% ------------------------------ PARAMETERS ------------------------------------
% Change these parameters to change which results to display
bitwidth      =  8; % Number of bits in input signal
fractionwidth =  4; % Number of fractional bits in input signal
output_bins   = 64; % Number of bins in the output
% ---------------------------- END PARAMETERS ----------------------------------

% Model parameters
N = 8;
M = output_bins;
r = 0.9990234375; % Same value as the Q12.10 fixed-point value

range = 0:output_bins-1;

% File parameters
% TODO: Insert path prefix
input_filename = sprintf('bwidth_%0d_fwidth_%0d_samples_%0d.csv', bitwidth, fractionwidth, N);
input_filepath = fullfile(input_file_prefix, input_filename);
input_datas = read_file(input_filepath);

output_filename_template = sprintf('bwidth_%0d_fwidth_%0d_bins_%0d_symbol_%s.csv', bitwidth, fractionwidth, output_bins, '%0d');

num_symbols = 4;
symbol_step = 2;
symbol_offset = 0;

figure;

for symbol_index = 0:num_symbols-1
    % Find the correct symbol to display
    symbol = symbol_index * symbol_step + symbol_offset;

    % Read the data from the SystemVerilog simulation
    output_filename = sprintf(output_filename_template, symbol);
    output_filepath = fullfile(output_file_prefix, output_filename);
    output_data = read_file(output_filepath);

    % Read the raw data from the input file
    input_data = input_datas(1, 1 + N*symbol:N + N*symbol);

    % Perform the rSDFT using the MATLAB model
    dft_data = rsdft(input_data, N, output_bins, r, range);

    % Find the magnitude for both models
    abs_data = abs(fftshift(output_data));
    abs_rsdft = abs(fftshift(dft_data));

    % Find the peak and its index
    [max_data, max_data_index] = max(abs_data);
    [max_rsdft, max_rsdft_index] = max(abs_rsdft);

    % Compute the difference between the spectra
    difference = abs(abs_data - abs_rsdft);

    % Plot the spectra and pinpoint the largest peak
    subplot(num_symbols, 2, 2*symbol_index + 1);
    hold on;
    plot(range, abs_data,  'Color', 'red', 'DisplayName', 'Implementation');
    plot(range, abs_rsdft, 'Color', 'blue', 'DisplayName', 'MATLAB model');
    if max_data_index == max_rsdft_index
        stem(max_data_index-1, max(max_data, max_rsdft), 'x', 'Color', 'green', 'HandleVisibility', 'off');
    else
        stem(max_data_index-1, max_data,  'x', 'Color', 'red', 'HandleVisibility', 'off');
        stem(max_rsdft_index-1, max_rsdft, 'x', 'Color', 'blue', 'HandleVisibility', 'off');
    end
    xlim([min(range), max(range)]);
    legend('location', 'best', 'Interpreter', 'latex');

    % Plot the difference between the magnitude spectra and pinpoint the largest peak
    subplot(num_symbols, 2, 2*symbol_index + 2);
    hold on;
    stem(range, difference, 'Color', 'black', 'DisplayName', '$|rsdft - model|$');
    if max_data_index == max_rsdft_index
        stem(max_data_index-1, difference(max_data_index), 'Color', 'green', 'HandleVisibility', 'off');
    else
        stem(max_data_index-1, difference(max_data_index), 'Color', 'red', 'HandleVisibility', 'off');
        stem(max_rsdft_index-1, difference(max_rsdft_index), 'Color', 'blue', 'HandleVisibility', 'off');
    end
    xlim([min(range), max(range)]);
    legend('location', 'best', 'Interpreter', 'latex');
end


% Read data from a file
function data = read_file(filepath)
    fid = fopen(filepath, 'r');

    if fid > 0
        header = strsplit(fgetl(fid), ' ');

        [lines, columns] = header{:};
        lines = str2double(lines);
        columns = str2double(columns);

        data = textscan(fid, '%f');

        data = data{:};

        if lines == size(data, 2) && columns == size(data, 1)
            data = transpose(data);
        end

        fclose(fid);
    else
        error('Could not open the file: %s', filepath);
    end
end
