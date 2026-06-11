function [txBits, txSym, txWave] = txChain(dataBits, preambleBits, p, txFilter)

txBits = [preambleBits(:); dataBits(:)];

txSym = qpskBitsToSymbols(txBits, p);

txWave = upfirdn(txSym, txFilter, p.sps, 1);

end