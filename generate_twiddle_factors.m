function twiddle_factors = generate_twiddle_factors(N, M)
% GENERATE_TWIDDLE_FACTORS  Generate M twiddle factors on the form W_M^Nk
%   twiddle_factors = GENERATE_TWIDDLE_FACTORS(1, M) generates the twiddle 
%   factors matching the expression exp( j * 2 * pi * k / M )
%
%   twiddle_factors = GENERATE_TWIDDLE_FACTORS(-N, M) generates the twiddle 
%   factors matching the expression exp( -j * 2 * pi * k * N / M )

% Create the vector of k (index) values
k = linspace(0, M-1, M);

twiddle_factors = exp(1j * 2 * pi * N * k / M);

end