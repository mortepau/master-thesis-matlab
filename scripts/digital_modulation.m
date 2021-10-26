%% ------- Constants -------

A_c = 1;
A_m = 1;

f_c = 200;
f_div = 100;

phi_c = 0;

samples = 10000;
symbols = 5;
samples_per_symbol = samples / symbols;
end_time = 0.1;

%% ------- Main -------

% Message signal
msg = binary_stream(A_m, symbols, samples_per_symbol);
% Carrier signal
carrier = to_signal(A_c, f_c, phi_c, samples, end_time);
carrier0 = to_signal(A_c, f_c-f_div, phi_c, samples, end_time);
carrier1 = to_signal(A_c, f_c+f_div, phi_c, samples, end_time);

bask = amplitude_shift_keying(msg, carrier);
bpsk = phase_shift_keying(msg, carrier);
bfsk = frequency_shift_keying(msg, carrier0, carrier1);

plot_triplet(msg, carrier, bask, end_time, true);
plot_triplet(msg, carrier, bpsk, end_time, true);
plot_triplet(msg, carrier, bfsk, end_time, true);

%% -------- Helper functions --------

function stream = binary_stream(A, num_samples, samples_per_symbol)
  stream = randi([0, 1], 1, num_samples);
  if num_samples == 10
    stream = [1, 1, 0, 1, 0, 1, 1, 0, 0, 1];
  elseif num_samples == 5
    stream = [1, 0, 1, 0, 1];
  end
  stream = A * repelem(stream, samples_per_symbol);
end

function signal = to_signal(A, f, phi, samples, end_time)
  w = 2 * pi * f;
  t = linspace(0, end_time, samples);

  signal = A * cos(w * t + phi);
end

function signal = amplitude_shift_keying(b, c)
  signal = b .* c(1:length(b));
end

function signal = phase_shift_keying(b, c)
  signal = c(1:length(b));
  for i = 1:length(b)
    if b(i) > 0
      signal(i) = -signal(i);
    end
  end
end

function signal = frequency_shift_keying(b, c0, c1)
  signal = zeros(1, length(b));
  for i = 1:length(b)
    if b(i) > 0
      signal(i) = c1(i);
    else
      signal(i) = c0(i);
    end
  end
end

function plot_signal(ax, end_time, signal, binary)
  x = linspace(0, end_time, length(signal));
  if binary
    stairs(ax, x, signal);
  else
    plot(ax, x, signal);
  end
  ylim([min(signal) - 0.1 * max(signal), 1.1 * max(signal)]);
  xlabel('time [t]');
  ylabel('Amplitude [V]');
end

function f = plot_triplet(msg, carrier, modulated, end_time, separate)
  if separate
    f = [];
    figure; f(1) = gca;
    plot_signal(f(1), end_time, msg, true);
    figure; f(2) = gca;
    plot_signal(f(2), end_time, carrier, false);
    figure; f(3) = gca;
    plot_signal(f(3), end_time, modulated, false);
  else
    f = figure;
    ax = subplot(2, 2, 1);
    plot_signal(ax, end_time, msg, true);
    ax = subplot(2, 2, 2);
    plot_signal(ax, end_time, carrier, false);
    ax = subplot(2, 2, [3, 4]);
    plot_signal(ax, end_time, modulated, false);
  end
end
