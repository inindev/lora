
clc
clear

test_inv = ~true;
if(test_inv)
    swap_iq = true;
    invert = ~swap_iq;
else
    swap_iq = false;
    invert = false;
end

% rx_sdr -g12 -f 910300000 -s 250000 -F CF32 /tmp/lora.raw
phy = LoraPhy(7, 125e3, '/tmp/lora.raw', 250e3, swap_iq);

[x, netid1, netid2] = phy.detect_preamble(1, 8, invert);
if(x < 1)
    fprintf('preamble not detected\n');
    return;
end
fprintf('netid1:%d  netid2:%d\n', netid1, netid2);

[pos_sfd, pos_hdr] = phy.detect_sfd(x, invert);
if(pos_sfd < 1)
    fprintf('sfd not detected\n');
    return;
end

[payload_len, cr, crc, is_valid] = phy.decode_header(pos_hdr, invert);
if(~is_valid)
    fprintf('header is invalid\n');
    return;
end
fprintf('header is valid - payload_len:%d  cr:%d  crc:%d\n', payload_len, cr, crc);

% payload symbols
payload_symbols = phy.decode_payload(pos_hdr, payload_len, invert);

% gray decode
payload_symbols_g = phy.gray_decode(payload_symbols, false);

% deinterleave
codewords = phy.diag_deinterleave(payload_symbols_g, phy.cr, phy.sf-2*phy.use_ldro);

% hamming decode
nibbles = phy.hamming_decode(codewords, phy.cr);
blen = uint8(floor(length(nibbles) / 2));
bytes = zeros(blen, 1);
for ii = 1:blen
    nn = ii*2;
    bytes(ii) = bitor(bitshift(nibbles(nn), 4), nibbles(nn-1));
end

% dewhiten
payload = phy.dewhiten(bytes(1:payload_len));

% display payload
fprintf('payload: ');
for ii = 1:length(payload)
    fprintf('%02x ', payload(ii))
end
fprintf('\n');
fprintf('csum1: %d\n', bytes(payload_len+1));
fprintf('csum2: %d\n', bytes(payload_len+2));

% phy.plot_symbols(1024, 13)
% phy.plot_symbol_pspec(3174);

