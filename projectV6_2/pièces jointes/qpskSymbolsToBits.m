function bits = qpskSymbolsToBits(sym, p)

demod = pskdemod(sym, p.M, p.phaseOffset, 'gray');

bitsMat = de2bi(demod, p.k, 'left-msb').';

bits = bitsMat(:);

end