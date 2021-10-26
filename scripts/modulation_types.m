% Debug, aka do not print
debug = true;
include_analog = false;
include_digital = true;

% Default to using latex interpretation of tick labels
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% Duration in seconds
duration = 0.1;

% Frequencies
f_s = 2000; % Sampling frequency
f_c = 200;  % Carrier frequency
f_m = 50;   % Message frequency

% Amplitudes
A_c = 2;
A_m = 2;

% Phases
P_c = 0;
P_m = pi / 4;

% Create the time vector
t = linspace(0, duration, f_s);

% Create the signals
signal_c = A_c * sin(2 * pi * f_c * t + P_c); % Carrier signal
signal_m = A_m * sin(2 * pi * f_m * t + P_m); % Message signal
signal_mI = A_m * sin(2 * pi * f_m * t) * cos(P_m); % Message signal, in-phase component
signal_mQ = A_m * cos(2 * pi * f_m * t) * sin(P_m); % Message signal, quadrature component

% Analog modulation

% - Frequency modulation
delta_f = f_c / 4;
signal_fm = A_c * cos( 2 * pi * f_c * t + P_c + A_m * delta_f / f_m * sin( 2 * pi * f_m * t + P_m) );

% - Amplitude modulation
modulation_sensitivity = 0.25;
signal_am = (1 + modulation_sensitivity * signal_m) .* signal_c;

% - Phase modulation
signal_pm = A_c * sin( 2 * pi * f_c * t + P_c + A_m * signal_m);

% Plot the signals

if include_analog

    % - Message signal
    figure;
    plot(t, signal_m);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_analog_msg_signal.png');
    end

    % - Message signal, IQ
    figure;
    hold on;
    plot(t, signal_mI);
    plot(t, signal_mQ);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_analog_msg_qi_signal.png');
    end

    % - Carrier signal
    figure;
    plot(t, signal_c);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_analog_carrier_signal.png');
    end

    % -- AM
    figure;
    plot(t, signal_am);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_analog_am_signal.png');
    end

    % -- FM
    figure;
    plot(t, signal_fm);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_analog_fm_signal.png');
    end

    % -- PM
    figure;
    plot(t, signal_pm);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_analog_pm_signal.png');
    end
    
end

% Digital modulation

% Variables
N_bits = 1; % Bits per symbol
N_symbol = 8; % Number of symbols
N = N_symbol * N_bits; % total number of bits
extension_factor = 40; % Number of bits to better interpolate the sine value

T_b = 1 / f_m;
T_s = T_b * N_bits;

E_s = 1 / N_bits * sum(0:2^N_bits-1);
E_b = E_s / N_bits;

t = linspace(0, N-1, N) * T_b;
t_extended = linspace(0, N-1, N * extension_factor) * T_b;

signal_mb = randi([0, 1], 1, N); % Message signal, binary
signal_mb_extended = repelem(signal_mb, extension_factor);

% - ASK modulation
prev_val = 0;
out = [];
for v = signal_mb
    if v
        if prev_val
            ts_start = 2*pi / (extension_factor - 1);
        else
            ts_start = 0;
        end
        ts = linspace(ts_start, 2 * pi, extension_factor);
        out = [out, cos(ts)]; % 2 * sqrt(E_b / T_b) * cos(ts)
    else
        out = [out, zeros(1, extension_factor)];
    end
    prev_val = v;
end
signal_ask = out;

% - FSK modulation
f_base = 3;
delta_f = 2;
f_1 = f_base - delta_f;
f_2 = f_base + delta_f;

prev_val = 0;
out = [];
ts_end = 0;
for v = signal_mb
    if v
        if ~ts_end || ~prev_val
            ts_end = 2 * pi * f_2 * ((extension_factor - 2) / (extension_factor - 1));
        end
        ts = linspace(0, ts_end, extension_factor);
    else
        if ~ts_end || prev_val
            ts_end = 2 * pi * f_1 * ((extension_factor - 2) / (extension_factor - 1));
        end
        ts = linspace(0, ts_end, extension_factor);
    end
    out = [out, cos(ts)]; % 2 * sqrt(E_b / T_b) * cos(ts)
    prev_val = v;
end
signal_fsk = out;

% - PSK modulation
prev_val = 0;
out = [];
for v = signal_mb
    if v
        phase = pi;
    else
        phase = 0;
    end
    ts_start = 2*pi / (extension_factor - 1) * (v == prev_val);
    ts = linspace(ts_start, 2 * pi, extension_factor);
    out = [out, cos(ts + phase)]; % 2 * sqrt(E_b / T_b) * cos(ts)
    prev_val = v;
end
signal_psk = out;

% - QAM modulation
signal_qam = signal_mb;

% Plot the signals

if include_digital
    
    xticks_ = t_extended(extension_factor * (1:N)) + (t_extended(2) - t_extended(1)) / 2;

    % -- Message signal
    figure;
    stairs(t_extended, signal_mb_extended);
    xticks(xticks_);
    %xticklabels(xlabels);
    ylim([-0.1, 1.1]);
    xlabel('Time [s]');
    ylabel('Binary value');
    yticks([0, 1]);
    grid on;
    ax = gca;
    ax.YGrid = 'off';
    if ~debug
        saveas(gfc, 'figures/back_digital_msg_signal.png');
    end

    % -- ASK
    figure;
    hold on;
    stairs(t_extended, signal_mb_extended);
    plot(t_extended, signal_ask);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    yticks([-1, 1])
    yticklabels({'-$2\sqrt{\frac{E_b}{T_b}}$', '$2\sqrt{\frac{E_b}{T_b}}$'});
    if ~debug
        saveas(gcf, 'figures/back_digital_ask_signal.png');
    end

    % -- FSK
    figure;
    hold on;
    stairs(t_extended, signal_mb_extended);
    plot(t_extended, signal_fsk);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_digital_fsk_signal.png');
    end

    % -- PSK
    figure;
    hold on;
    stairs(t_extended, signal_mb_extended);
    plot(t_extended, signal_psk);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    if ~debug
        saveas(gcf, 'figures/back_digital_psk_signal.png');
    end

    % -- QAM
    if ~debug
    figure;
    stairs(t, signal_qam);
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    
        saveas(gcf, 'figures/back_digital_qam_signal.png');
    end

end