function compute_rsdft(bitwidth, fractionwidth, output_bins, ebn, packet)
% Compute the rSDFT of the values in the file specified by the parameters

% ------------------------------ PARAMETERS ------------------------------------
% Change these parameters to change which results to display
% bitwidth      =  8; % Number of bits in input signal
% fractionwidth =  4; % Number of fractional bits in input signal
% output_bins   = 64; % Number of bins in the output
% ebn           = 14; % EbN0 value to use
% packet        =  1; % Packet to use
% ---------------------------- END PARAMETERS ----------------------------------

if nargin < 5
  error('Not enough input arguments');
end

% Fixed-Point parameters
fi_type = numerictype( ...
  'DataType',       'Fixed', ...
  'DataTypeMode',   'Fixed-point: binary point scaling', ...
  'Signedness',     'Signed', ...
  'WordLength',     bitwidth, ...
  'FractionLength', fractionwidth ...
);
fi_math = fimath( ...
  'RoundingMethod',        'Zero', ...
  'OverflowAction',        'Saturate', ...
  'ProductMode',           'SpecifyPrecision', ...
  'ProductWordLength',     fi_type.WordLength, ...
  'ProductFractionLength', fi_type.FractionLength, ...
  'SumMode',               'SpecifyPrecision', ...
  'SumWordLength',         fi_type.WordLength, ...
  'SumFractionLength',     fi_type.FractionLength, ...
  'CastBeforeSum',         true ...
);

% Model parameters
N = 8;
M = output_bins;
r = 0.9990234375; % Same value as the Q12.10 fixed-point value

% File parameters
script_path = fileparts(mfilename('fullpath'));
input_file_prefix = '../data/receiver/';
input_filename = sprintf('bwidth_%0d_fwidth_%d_N_%d_ebn_%d_packet_%d.csv', bitwidth, fractionwidth, N, ebn, packet);
input_filepath = fullfile(script_path, input_file_prefix, input_filename);
input_datas = read_file(input_filepath);

output_file_prefix = '../data/rsdft/';
output_filename_template = sprintf('bwidth_%d_fwidth_%d_N_%d_bins_%d_ebn_%d_packet_%d_symbol_%s.csv', bitwidth, fractionwidth, N, output_bins, ebn, packet, '%0d');


fprintf("Computing the %d-point rSDFT of the data defined by the parameters:\n", M);
fprintf('Qi.f  : Q%d.%d\n', bitwidth, fractionwidth);
fprintf('N     : %d\n', N);
fprintf('M     : %d\n', M);
fprintf('EbN0  : %d\n', ebn);
fprintf('Packet: %d\n', packet);

num_indexes = length(input_datas) / N - 1;

for index = 0:num_indexes
    % Samples for a single symbol
    data = input_datas(1 + index*N : N + index*N);

    % Convert it to Fixed-Point
    fi_data = fi(data, fi_type, fi_math);

    % Compute the rSDFT values
    rsdft_data = rsdft(fi_data, N, M, r, 0:M-1);

    % Open the file and write the content of data to it
    output_filename = sprintf(output_filename_template, index);
    output_filepath = fullfile(script_path, output_file_prefix, output_filename);
    write_file(output_filepath, rsdft_data);

    if index == 0
      fprintf(1, 'Progress: %4d / %-4d ( %3d%% )', index, num_indexes, index / num_indexes * 100);
    else
      fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%4d / %-4d ( %3d%% )', index, num_indexes, floor(index / num_indexes * 100));
    end
end
fprintf('\n--- FINISHED ---\n');

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

% Write data to a file
function write_file(filepath, data)
    fid = fopen(filepath, 'w+');

    spec = "%.16f+%.16fi\n";

    if fid > 0
        % Print the size of the data
        fprintf(fid, '%d %d\n', size(data));

        % Print the actual data
        if isreal(data)
            fprintf(fid, spec, data');
        else
            fprintf(fid, spec, [real(data(:)), imag(data(:))].');
        end

        fclose(fid);
    else
        error('Could not open the file: %s', filepath);
    end
end

