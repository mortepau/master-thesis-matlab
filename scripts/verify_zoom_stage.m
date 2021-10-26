% Verify the functionality of the implementation of the zoom stage
% This is done by comparing the generated BoiFinal from the hardware
% implementation with the simulated MATLAB model.

% ------------------------------ PARAMETERS ------------------------------------
% Change these parameters to change which results to display
bitwidth         =  8; % Number of bits in input signal
fractionwidth    =  4; % Number of fractional bits in input signal
boi_final_length = 16; % Number of indexes in the BoI_Final
% ---------------------------- END PARAMETERS ----------------------------------

% Model parameters
N = 8;
I = 8;
r = 0.9990234375; % Same value as the Q12.10 fixed-point value
ebn = 14;
packet = 2;

% File parameters

input_filename = sprintf('bwidth_%0d_fwidth_%0d_symbol_%0d_ebn_%0d_packet_%0d.csv', bitwidth, fractionwidth, N, ebn, packet);
input_filepath = fullfile(input_file_prefix, input_filename);
input_datas = read_file(input_filepath, '%f');

% Read the data from the SystemVerilog simulation
output_filename = sprintf('symbol_%0d_interpolation_%0d_boi_%0d_ebn_%0d_packet_%0d.csv', N, I, boi_final_length, ebn, packet);
output_filepath = fullfile(output_file_prefix, output_filename);
output_data = read_file(output_filepath, '%d');

% Read the raw data from the input file
% input_data = input_datas(1, 1 + N*symbol:N + N*symbol);
input_data = input_datas;

% Find the BoI from the MATLAB model
boi_final = zoom_stage(input_datas, N, I, boi_final_length, false);

% Compute the rSDFT spectrum
rsdft_data = abs(rsdft(input_data, N, N*I, r, 0:N*I-1));

% Find the overlapping indexes
[~, overlapping_indexes] = intersect(boi_final, output_data+1);
overlapping_indexes = sort(overlapping_indexes);
overlapping = boi_final(overlapping_indexes);
nonoverlapping_input = setdiff(boi_final, output_data+1);
nonoverlapping_output = setdiff(output_data+1, boi_final);

figure;

hold on;
plot(0:N*I-1, rsdft_data, 'Color', '#000', 'DisplayName', '$X_n(k)$');
plot(overlapping-1, rsdft_data(overlapping), 'o', 'Color', '#F0F', 'MarkerFaceColor', '#F0F', 'DisplayName', "overlapping k's");
plot(nonoverlapping_input-1, rsdft_data(nonoverlapping_input), 'o', 'Color', '#F00', 'MarkerFaceColor', '#F00', 'DisplayName', "non-overlapping MATLAB k's");
plot(nonoverlapping_output-1, rsdft_data(nonoverlapping_output), 'o', 'Color', '#00F', 'MarkerFaceColor', '#00F', 'DisplayName', "non-overlapping SystemVerilog k's");

xlim([0, N*I-1]);
legend('location', 'best', 'Interpreter', 'latex');

% Read data from a file
function data = read_file(filepath, format)
    fid = fopen(filepath, 'r');

    if fid > 0
        header = strsplit(fgetl(fid), ' ');

        [lines, columns] = header{:};
        lines = str2double(lines);
        columns = str2double(columns);

        data = textscan(fid, format);

        data = data{:};

        if lines == size(data, 2) && columns == size(data, 1)
            data = transpose(data);
        end

        fclose(fid);
    else
        error('Could not open the file: %s', filepath);
    end
end
