function dft = goertzel(samples, k)

% Make it 0-indexed
k = k - 1;

N = length(samples);

% Extend the samples with a zero at the tail
samples = [ 0, samples, 0 ];

w0 = 2 * pi * k / N;

expw0 = exp( -j * w0 );
cosw0 = 2 * cos( w0 );

expk = exp( -j * 2 * pi * k );

ykn = zeros(3, length(k));

for i = 2:N+2
  ykn(3, :) = samples(i) - expw0 * samples(i-1) + cosw0 .* ykn(2, :) - ykn(1, :);
  ykn = [ ykn(2:3, :); zeros(1, length(k)) ];
end
dft = ykn(2,:);

end
