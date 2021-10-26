% Create a plot showcasing the effect the r-factor in the rSDFT has on the
% spectra

N = 8;
I = 8;

fc = 400;
dt = 1 / fc;
stopTime = 0.02;
t = 0:dt:stopTime-dt;

f = 105;
signal = cos(2 * pi * f * t) + 1j * sin(2 * pi * f * t);

figure;
hold on;
grid on;
grid minor;

plot(1:N*I, abs(fft(signal, N*I)), 'LineWidth', 1, 'MarkerFaceColor', 'red', 'DisplayName', '64-point FFT');
for r = [1, 0.9, 0.8, 0.7]
  res = rsdft(signal, N, N*I, r, 0:N*I-1);
  plot(1:N*I, abs(res), 'LineWidth', 1, 'DisplayName', sprintf('64-point rSDFT r = %0.2f', r));
end

xlabel('f [Hz]')
ylabel('S(f)')
xlim([1, N*I])
xticks(1:(N-1):N*I)
labels_as_num = num2cell(round(linspace(0, fc, 10)));
xticklabels(labels_as_num)
legend();
