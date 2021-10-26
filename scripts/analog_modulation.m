%% ------- Constants -------

A_c = 3;
A_m = 5;

f_c = 200;
f_m =  15;
f_dev = 33;

phi_c = 0;
phi_m = 0;

samples = 10000;
end_time = 0.1;

%% ------- Main -------

% Message signal
msg = to_signal(A_m, f_m, phi_m, samples, end_time);
msg_d = msg > 0;
% Carrier signal
carrier = to_signal(A_c, f_c, phi_c, samples, end_time);

% Amplitude modulated signal
am = amplitude_modulation(A_m, f_m, phi_m, A_c, f_c, phi_c, samples, end_time);
% Phase modulated signal
pm = phase_modulation(A_m, f_m, phi_m, A_c, f_c, phi_c, f_dev, samples, end_time);
% Frequency modulated signal
fm = frequency_modulation(A_m, f_m, phi_m, A_c, f_c, phi_c, f_dev, samples, end_time);

plot_triplet(msg, carrier, am, end_time, true);
figure; ax = gca;
plot_signal(ax, end_time, msg_d);
figure; ax = gca;
plot_signal(ax, end_time, pm);
figure; ax = gca;
plot_signal(ax, end_time, fm);


%% -------- Helper functions --------

function signal = to_signal(A, f, phi, samples, end_time)
  w = 2 * pi * f;
  t = linspace(0, end_time, samples);

  signal = A * cos(w * t + phi);
end

function signal = amplitude_modulation(A_m, f_m, phi_m, A_c, f_c, phi_c, samples, end_time)
  w_m = 2 * pi * f_m;
  w_c = 2 * pi * f_c;
  t = linspace(0, end_time, samples);
  A = A_c * ( 1 + A_m / A_c ) * cos( w_m * t + phi_m);
  signal = A .* cos( w_c * t + phi_c);
end

function signal = phase_modulation(A_m, f_m, phi_m, A_c, f_c, phi_c, f_dev, samples, end_time)
  w_m = 2 * pi * f_m;
  w_c = 2 * pi * f_c;
  t = linspace(0, end_time, samples);
  signal = A_c * cos( w_c * t + f_dev * A_m * cos( w_m * t + phi_m ) );
end

function signal = frequency_modulation(A_m, f_m, phi_m, A_c, f_c, phi_c, f_dev, samples, end_time)
  w_m = 2 * pi * f_m;
  w_c = 2 * pi * f_c;
  beta = f_dev * A_m / f_m;
  t = linspace(0, end_time, samples);
  signal = A_c * cos( w_c * t + beta * sin( w_m * t + phi_m ) );
end

function plot_signal(ax, end_time, signal)
  x = linspace(0, end_time, length(signal));
  plot(ax, x, signal);
  ylim([min(signal) - 0.1 * max(signal), 1.1 * max(signal)]);
  xlabel('time [t]');
  ylabel('Amplitude [V]');
end

function f = plot_triplet(msg, carrier, modulated, end_time, separate)
  if separate
    f = [];
    figure; f(1) = gca;
    plot_signal(f(1), end_time, msg);
    figure; f(2) = gca;
    plot_signal(f(2), end_time, carrier);
    figure; f(3) = gca;
    plot_signal(f(3), end_time, modulated);
  else
    f = figure;
    ax = subplot(2, 2, 1);
    plot_signal(ax, end_time, msg);
    ax = subplot(2, 2, 2);
    plot_signal(ax, end_time, carrier);
    ax = subplot(2, 2, [3, 4]);
    plot_signal(ax, end_time, modulated);
  end
end
