function packet = generate_packet(packet_length, preamble_length)
% GENERATE_PACKET  Create a packet to be used for the demodulator.
%   data = GENERATE_PACKET(N, M) generates a 1xN vector where the M first
%   values are an alternating sequence of 1's and 0's, and the N-M-1 last
%   values is a pseudorandom sequence of 1's and 0's.
%
%   GENERATE_PACKET will raise an error if preamble_length is greater than
%   packet_length.

    % Assert that the preamble_length is not longer than the packet_length.
    if preamble_length > packet_length
        error('preamble_length cannot be greater than packet_length')
    end

    % Create a vector of alternating 1's and 0's
    preamble = [ ones(1, floor(preamble_length/2)); zeros(1, floor(preamble_length/2)) ];
    preamble = preamble(:)';

    % If the preamble_length is odd then we append a 1 at the end of it, as
    % the created alternating series ends with a 0.
    if mod(preamble_length, 2)
        preamble(end+1) = 1;
    end

    message = randi([0, 1], 1, packet_length - preamble_length - 1);

    % Assumes that the delimiter between the preamble and the message is a
    % bit such that the alternating pattern is broken
    delimiter = preamble(end);

    packet = [ preamble, delimiter, message ];
end
