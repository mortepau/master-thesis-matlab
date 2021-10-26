function dft = rsdft(samples, N, M, r, k)
% RSDFT Calculate the DFT values using the rSDFT algorithm
%   dft = RSDFT(samples, N, M, r, num_shifts, k) calculates the DFT values
%   for samples using the rSDFT algorithm with M twiddle factors for a
%   window of length N over num_shifts for the bins in k (0-indexed).

W_k = generate_twiddle_factors(1, M);
W_Nk = generate_twiddle_factors(-N, M);

dft = zeros(1, length(k));

% Append N samples at the start so we "remove" 0 for the N first samples
samples = [zeros(1, N), samples];

for shift = 1:length(samples)-N
    dft = W_k(k+1) .* ( r * dft - r^N * samples(shift) + W_Nk(k+1) .* samples(shift + N) );
end

end
