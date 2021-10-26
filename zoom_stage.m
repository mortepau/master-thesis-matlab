function BoI_final = zoom_stage(samples, N, I, BoI_final_length, include_graphics)
% ZOOM_STAGE Perform the zoom stage in the demodulator algorithm
% BoI_final = ZOOM_STAGE(samples, N, I, BoI_final_length, false) Finds the
% Bins of Interest for the N*I-point DFT of samples and outputs BoI_final
% with a length of BoI_final_length.

% TODO: Check if the rSDFT should only calculate on BoI or all bins in the DFT

BoI = 1:N;

if include_graphics
    figure;
end

% Iterate over all the gamma values
for gamma = 0:log2(I)

    dft = zeros(2 * N, N * 2^(gamma));

    % Calculate the DFT of the window using the BoI
    parfor i = 1:2*N
        window = samples(i : N - 1 + i);

        % Calculate the DFT of the window
        dft(i, :) = abs(rsdft(window, N, N * 2^(gamma), 1, 0:N * 2^(gamma)-1));
    end

    % Sum the bins for all k's across i's
    magnitude = sum(dft);

    % Find the magnitudes for the bins in BoI
    magnitude_BoI = magnitude(BoI);

    if include_graphics
        subplot(log2(I)+2, 1, gamma+1);
        plot(0 : N * 2 ^ gamma - 1, magnitude);
        hold on;
        plot(BoI-1, magnitude_BoI, 'o');
        xlim([0, N * 2 ^ gamma - 1]);
        xlabel('k');
        ylabel('$|X(f)|$', 'Interpreter', 'latex');
    end

    % Find the index for the value with the absolute maximal value
    [~, k_c_index] = max(magnitude_BoI);
    k_c = BoI(k_c_index);

    % Estimate the bin for the carrier frequency in the next stage
    % Subtract 1 from k_c to have it 0-indexed instead of 1-indexed
    k_c_estimate = 2 * (k_c - 1);

    % These are for a 0-indexed array
    a = mod(k_c_estimate - I / 2, N * 2^(gamma + 1));
    b = mod(k_c_estimate + I / 2 - 1, N * 2^(gamma + 1));

    % Create the indexes matching the frequency bin, and the "one"-bin and "zero"-bin
    if a < b
        BoI = a:b;
    else
        % Wrap around so that the bin order is preserved
        BoI = [a:N * 2^(gamma+1) - 1, 0:b];
    end

    % Make it 1-indexed
    BoI = BoI + 1;

end

% Same process as for a and b, but use BoI_final_length instead of I
d = mod((k_c - 1) - BoI_final_length / 2, N * I);
e = mod((k_c - 1) + BoI_final_length / 2 - 1, N * I);

if d < e
	BoI_final = d:e;
else
	BoI_final = [d:N * I - 1, 0:e];
end
BoI_final = BoI_final + 1;

if include_graphics
    subplot(log2(I)+2, 1, log2(I)+2);
    plot(0 : N * I - 1, magnitude);
    hold on;
    plot(BoI_final-1, magnitude(BoI_final), 'o');
    xlim([0, N * I - 1]);
    xlabel('k');
    ylabel('$|X(f)|$', 'Interpreter', 'latex');
end

end
