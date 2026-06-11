function sym = qpskBitsToSymbols(bits, p)

bits = bits(:);

if mod(length(bits), p.k) ~= 0
    error('Number of bits must be multiple of p.k.');
end

intData = bi2de(reshape(bits, p.k, []).', 'left-msb');

sym = pskmod(intData, p.M, p.phaseOffset, 'gray');
sym = sym(:);

end