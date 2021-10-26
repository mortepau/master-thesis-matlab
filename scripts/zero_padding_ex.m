% Create a plot showcasing the effect zero-padding has on the spectra,
% and where the original samples will be in the zero-padded spectra

N = 8;
I = 8;

fc = 400;
dt = 1 / fc;
stopTime = 0.02;
t = 0:dt:stopTime-dt;

f = 177;
signal = cos(2 * pi * f * t) + 1j * sin(2 * pi * f * t);

zero_padding = zeros(1, (I-1)*N);

non_zero_padded_fft  = fft(signal(1:N), N);
zero_padded_fft = fft([signal(1:N), zero_padding], N*I);

nzp_spectrum = abs(non_zero_padded_fft);
zp_spectrum = abs(zero_padded_fft);

figure;
hold on;
grid on;
grid minor;

plot(1:N, nzp_spectrum, '-o', 'LineWidth', 1, 'MarkerEdgeColor', 'blue');
plot(1:N*I, zp_spectrum, 'LineWidth', 1, 'MarkerFaceColor', 'red');
plot(1:I:N*I, nzp_spectrum, 'd', 'LineWidth', 1, 'MarkerEdgeColor', 'magenta');
for x = [4, 5, 6]
quiver(x, nzp_spectrum(x), (x-1)*(N - 1) - 0.5, 0, 0, 'Color', 'black', 'MaxHeadSize', 0.04)
end

xlabel('f [Hz]')
ylabel('S(f)')
xlim([1, N*I])
xticks(1:(N-1):N*I)
labels_as_num = num2cell(round(linspace(0, fc, 10)));
xticklabels(labels_as_num)
legend({ 'original signal', 'zero-padded signal', 'original signal, spaced' }, 'Location', 'best')
