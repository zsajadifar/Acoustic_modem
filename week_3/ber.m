function BER = ber(x , x_demod)

BER = sum(abs(x_demod-x));
BER = BER / length(x);