function [ k_0, k_1, timing ] = window_alignment_stage(samples, N, I, BoI, synchronization_length, num_delays, include_graphics)

if include_graphics
    figure;
end

R    = zeros(1, num_delays);
k_Es = zeros(1, num_delays);
k_Os = zeros(1, num_delays);

dft_even = zeros(synchronization_length / 2, length(BoI));
dft_odd  = zeros(synchronization_length / 2, length(BoI));

zero_padding = zeros(1, (I - 1) * N);

for delay = 1:num_delays
    parfor sync = 0:synchronization_length / 2 - 1
        samples_even = samples( (2 * sync)     * N + delay : (2 * sync)     * N + delay + N - 1);
        samples_odd  = samples( (2 * sync + 1) * N + delay : (2 * sync + 1) * N + delay + N - 1);

        dft_even(sync + 1, :) = abs(goertzel([ samples_even, zero_padding ], BoI));
        dft_odd(sync + 1, :)  = abs(goertzel([ samples_odd,  zero_padding ], BoI));
    end

    sum_even = sum(dft_even);
    sum_odd = sum(dft_odd);

    [~, k_E] = max(sum_even);
    [~, k_O] = max(sum_odd);

    R(delay) = ( sum_even(k_E) - sum_even(k_O) ) + ( sum_odd(k_O) - sum_odd(k_E) );
    k_Es(delay) = BoI(k_E);
    k_Os(delay) = BoI(k_O);
end

[~, timing] = max(R);

even_in_upper_half = k_Es(timing) > N * I / 2;
odd_in_upper_half = k_Os(timing) > N * I / 2;

if (even_in_upper_half && odd_in_upper_half) || (~even_in_upper_half && ~odd_in_upper_half)
    k_0 = max(k_Es(timing), k_Os(timing));
    k_1 = min(k_Es(timing), k_Os(timing));
else
    k_0 = min(k_Es(timing), k_Os(timing));
    k_1 = max(k_Es(timing), k_Os(timing));
end

end
