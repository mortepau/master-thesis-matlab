function [BER, BER_final] = tb_demodulator(use_fixed_point, use_amplitude, bitwidth, fraction_width, ebns, num_packets)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          Setup                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Debugging
debug_zoom_stage             = false;
debug_window_alignment_stage = false;

% Inject inputs at any stage in the pipeline
inject_at_receiver               = false;
inject_at_zoom_stage             = false;
inject_at_window_alignment_stage = false;
inject_at_data_detection_phase   = false;
inject_capture = inject_at_receiver || inject_at_zoom_stage || inject_at_window_alignment_stage;

% Eject outputs at any stage in the pipeline
eject_at_transmitter            = false;
eject_at_receiver               = false;
eject_at_zoom_stage             = false;
eject_at_window_alignment_stage = false;
eject_at_data_detection_phase   = false;
eject_at_error_analysis         = false;
eject_capture = eject_at_transmitter || eject_at_receiver || eject_at_window_alignment_stage;

transmitter_template            = '%sN_%d_ebn_%d_packet_%d.csv';
receiver_template               = '%sN_%d_ebn_%d_packet_%d.csv';
zoom_stage_template             = '%sN_%d_I_%d_boi_%d_ebn_%d_packet_%d.csv';
window_alignment_stage_template = '%sN_%d_I_%d_boi_%d_L_%d_sync_%d_ebn_%d_packet_%d.csv';
data_detection_phase_template   = '%sN_%d_I_%d_boi_%d_L_%d_sync_%d_P_%d_ebn_%d_packet_%d.csv';
error_analysis_template         = '%sN_%d_I_%d_boi_%d_L_%d_sync_%d_P_%d_ebn_%d_packet_%d.csv';

% Fixed-Point parameters
fi_type = numerictype( ...
  'DataType',       'Fixed', ...
  'DataTypeMode',   'Fixed-point: binary point scaling', ...
  'Signedness',     'Signed', ...
  'WordLength',     bitwidth, ...
  'FractionLength', fraction_width ...
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


% System parameters
symbol_length          =  8;            % Samples per symbol
synchronization_length = 16;            % Number of symbols used in synchronization
BoI_length             = 16;            % Length of Bins of Interest ( BoI )
delay                  =  0;            % Value in range (0, N-1 = symbol_length - 1 )
num_window_delays      = symbol_length; % Number of delay values tested in synchronization phase
interpolation_factor   = 8;             % Interpolation factor

% Signal parameters
r_b       = 100;                                     % Bit rate
T_b       = 1 / r_b;                                 % Bit period
F_s       = symbol_length * r_b;                     % Sampling frequency
T_s       = 1 / F_s;                                 % Sampling period
f_dev     = r_b;                                     % FSK frequency separation
f_offset  = 35;                                      % Frequency offset value
amplitude = 2 ^ (bitwidth - fraction_width - 7);
% amplitude = 2 ^ (bitwidth - fraction_width - 1) - 1 % Amplitude of received signal

% Packet data parameters
% num_packets     =   2; % Number of packets to calcuate BER from
packet_length   = 100; % Symbols in the packet
preamble_length =   32; % Symbols in the preamble

% Noise parameters
ebn_start =  1;
ebn_end   = 14;
% ebns      = ebn_start:ebn_end; % E_b / N_0 values

% Error Analysis
analysis_start_index = preamble_length + 1;
analysis_end_index   = packet_length - 1;
analysis_range       = analysis_end_index - analysis_start_index + 1;

% Matrix holding all the calculated BER values
BER = zeros(length(ebns), num_packets);

% Transmitter objects
modulator = comm.FSKModulator( ...
    'ModulationOrder',     2, ...                            % Number of frequencies in signal
    'BitInput',            true, ...                         % The input are bits not integers
    'FrequencySeparation', f_dev, ...                        % Separation between successive tones
    'SamplesPerSymbol',    symbol_length ...                 % Number of samples per symbol
);

% Channel objects
channel_offset = comm.PhaseFrequencyOffset( ...
    'FrequencyOffsetSource', 'Property', ...
    'FrequencyOffset', f_offset, ...                         % Frequency offset which is applied
    'SampleRate',      F_s ...                               % Sample rate of input signal
);
channel_noise  = comm.AWGNChannel( ...
    'NoiseMethod',      'Signal to noise ratio (Eb/No)', ... % Which metric is used for the noise
    'SignalPower',      1, ...                               % Signal power in Watts ( 1W = 0dBW )
    'SamplesPerSymbol', symbol_length ...                    % Number of samples per symbol
);

% Notation matching paper
N = symbol_length;
I = interpolation_factor;
L = preamble_length;

% Check if the Parallel Computing Toolbox can be used
if license('test', 'distrib_computing_toolbox')
  pool = gcp('nocreate');
  if isempty(pool)
      pool = parpool('local');
  end
end

% Suffix for filename when using Fixed-Point
file_prefix = '';
if use_fixed_point
  if use_amplitude
    file_prefix = sprintf('bwidth_%d_fwidth_%d_amp_1_', bitwidth, fraction_width);
  else
    file_prefix = sprintf('bwidth_%d_fwidth_%d_amp_0_', bitwidth, fraction_width);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        Program                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = datetime('now', 'Format', 'hh:mm:ss.SSS');

fprintf('[%s] - Starting execution\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));

for ebn_index = 1:length(ebns)

    start_time_ebn = datetime('now', 'Format', 'ss.SSS');

    for packet_index = 1:num_packets
        packet_tic = tic;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                  Transmitter                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Create a new vector of data representing a packet
        if inject_capture
          packet_at_transmitter = inject('transmitter', sprintf(transmitter_template, file_prefix, N, ebns(ebn_index), packet_index));
        else
          packet_at_transmitter = generate_packet(packet_length, L);
        end

        % Transmit the packet
        packet_after_modulation = modulator(packet_at_transmitter');

        if eject_at_transmitter || eject_capture
            eject('transmitter', packet_at_transmitter, sprintf(transmitter_template, file_prefix, N, ebns(ebn_index), packet_index));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                    Channel                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Update the E_b / N_0 of the AWGNChannel object
        channel_noise.EbNo = ebns(ebn_index);

        % Modify the packet to simulate a channel
        packet_in_channel  = channel_offset(packet_after_modulation);
        packet_at_receiver = channel_noise(packet_in_channel)';

        % Multiply by a factor to simulate the amplitude
        % Note: This is done here to apply amplitude to the noise as well
        if use_amplitude
          packet_at_receiver = amplitude * packet_at_receiver;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                   Receiver                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Need to inject a "known" packet if we want to reuse results later in the pipeline
        if inject_at_receiver || inject_at_zoom_stage || inject_at_window_alignment_stage || inject_at_data_detection_phase
            packet = inject('receiver', sprintf(receiver_template, file_prefix, N, ebns(ebn_index), packet_index));

            % Convert the values to Fixed-Point
            if use_fixed_point && license('test', 'fixed_point_toolbox')
                packet = fi(packet, fi_type, fi_math);
            end
        else
            % Delay the packet, by removing the delay first samples and append delay 0's at the tail
            packet = [ packet_at_receiver(1, delay + 1 : end), zeros(1, delay) ];
            % fprintf('[%s] - Packet received\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));

            % Convert the values to Fixed-Point
            if use_fixed_point && license('test', 'fixed_point_toolbox')
                packet = fi(packet, fi_type, fi_math);
            end

            if eject_at_receiver || eject_capture
                eject('receiver', packet, sprintf(receiver_template, file_prefix, N, ebns(ebn_index), packet_index));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                  Zoom Stage                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~inject_at_window_alignment_stage
          % fprintf('[%s] - Entering zoom stage\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));
          % zoom_stage_tic = tic;
          if inject_at_zoom_stage
              BoI_final = inject('zoom_stage', sprintf(zoom_stage_template, file_prefix, N, I, BoI_length, ebns(ebn_index), packet_index));
          else
              BoI_final = zoom_stage(packet, N, I, BoI_length, debug_zoom_stage);

              if eject_at_zoom_stage
                  eject('zoom_stage', BoI_final, sprintf(zoom_stage_template, file_prefix, N, I, BoI_length, ebns(ebn_index), packet_index));
              end
          end
          % fprintf('[%s] - Leaving zoom stage\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));
          % fprintf('[%s] - Zoom stage completed in %gs\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), toc(zoom_stage_tic));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                            Window Alignment Stage                            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % fprintf('[%s] - Entering window alignment stage\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));
        % window_alignment_tic = tic;
        if inject_at_window_alignment_stage
            [k_0, k_1, timing] = inject('window_alignment_stage', sprintf(window_alignment_stage_template, file_prefix, N, I, BoI_length, L, synchronization_length, ebns(ebn_index), packet_index));
        else
            [k_0, k_1, timing] = window_alignment_stage(packet, N, I, BoI_final, synchronization_length, num_window_delays, debug_window_alignment_stage);

            if eject_at_window_alignment_stage
                eject('window_alignment_stage', [k_0, k_1, timing], sprintf(window_alignment_stage_template, file_prefix, N, I, BoI_length, L, synchronization_length, ebns(ebn_index), packet_index));
            end
        end
        % fprintf('[%s] - Leaving window alignment stage\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));
        % fprintf('[%s] - Window alignment stage completed in %gs\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), toc(window_alignment_tic));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                Data Detection                                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % fprintf('[%s] - Entering data detection stage\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));
        data_detection_tic = tic;
        if inject_at_data_detection_phase
            estimated_packet = inject('data_detection_phase', sprintf(data_detection_phase_template, file_prefix, N, I, BoI_length, L, synchronization_length, num_packets, ebns(ebn_index), packet_index));
        else
            estimated_packet = zeros(1, packet_length);
            goertzel_values = zeros(packet_length, 2);

            % Using a parfor here to decrease the runtime as everything can be
            % done in parallel
            for index = 0:packet_length - preamble_length - 1 - 1
                symbol_tic = tic;

                % Extract the symbol we want to detect value for
                symbol = packet( N * ( L + index ) + timing : N * ( L + index ) + timing + N - 1);

                % Add zero-padding so the total length is N * I
                zero_padding = zeros(1, (I - 1) * N);

                % Calculate the magnitude of the 0-bin and 1-bin
                magnitudes = abs(goertzel([ symbol, zero_padding ], [k_0, k_1]));

                % Split the vector into its respective component
                zero = magnitudes(1);
                one = magnitudes(2);

                % Estimate the value of the symbol
                estimated_packet(index + 1) = one > zero;

                % fprintf('[%s] - Symbol estimation completed in %gs\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), toc(symbol_tic));
            end

            if eject_at_data_detection_phase
                eject('data_detection_phase', estimated_packet, sprintf(data_detection_phase_template, file_prefix, N, I, L, synchronization_length, num_packets, ebns(ebn_index), packet_index));
            end
        end
        % fprintf('[%s] - Leaving data detection stage\n', datetime('now', 'Format', 'HH:mm:ss.SSS'));
        % fprintf('[%s] - Data detection stage completed in %gs\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), toc(data_detection_tic));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                Error Analysis                                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Calculate different kinds of error to compensate for the possible
        % misalignment of detected data by one or two symbols.
        % This would in practice be done by a frame synchronizer, but this
        % is out of the scope for this testbench

        error_no_shift           = sum(packet_at_transmitter(analysis_start_index + 0 : analysis_end_index + 0) ~= estimated_packet(1 + 0 : analysis_range + 0));
        error_single_left_shift  = sum(packet_at_transmitter(analysis_start_index + 1 : analysis_end_index + 1) ~= estimated_packet(1 + 0 : analysis_range + 0));
        error_single_right_shift = sum(packet_at_transmitter(analysis_start_index + 0 : analysis_end_index + 0) ~= estimated_packet(1 + 1 : analysis_range + 1));
        error_double_right_shift = sum(packet_at_transmitter(analysis_start_index + 0 : analysis_end_index + 0) ~= estimated_packet(1 + 2 : analysis_range + 2));

        BER(ebn_index, packet_index) = min([ error_no_shift, error_single_left_shift, error_single_right_shift, error_double_right_shift ]) / analysis_range;

        if eject_at_error_analysis
            eject('error_analysis', BER(ebn_index,:), sprintf(error_analysis_template, file_prefix, N, I, BoI_length, L, synchronization_length, num_packets, ebns(ebn_index), packet_index));
        end

        % fprintf('[%s] - Packet completed in %gs\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), toc(packet_tic));
        % Progress bar
        pb_current = (ebn_index-1) * num_packets + packet_index;
        pb_total = length(ebns) * num_packets;
        if ebn_index == 1 && packet_index == 1
          fprintf(1, 'Progress: %4d / %-4d ( %3d%% )', pb_current, pb_total, floor(pb_current / pb_total * 100));
        else
          fprintf(1, '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%4d / %-4d ( %3d%% )', pb_current, pb_total, floor(pb_current / pb_total * 100));
          if ebn_index == length(ebns) && packet_index == num_packets
            fprintf(1, '\n');
          end
        end
    end

    end_time_ebn = datetime('now', 'Format', 'ss.SSS');

    % fprintf('[%s] - E_b/N_0 execution time [%02d/%02d]: %s (Total: %s)\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), ebn_index, length(ebns), duration(end_time_ebn - start_time_ebn, 'Format', 'mm:ss.SSS'), duration(end_time_ebn - start_time, 'Format', 'mm:ss.SSS'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      Finalizing                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum the BER values for each E_b / N_0 value
BER_final = sum(BER, 2) / num_packets;

end_time = datetime('now', 'Format', 'hh:mm:ss.SSS');
fprintf('[%s] - Total execution time: %s\n', datetime('now', 'Format', 'HH:mm:ss.SSS'), duration(end_time - start_time, 'Format', 'mm:ss.SSS'));
end
