%
% lora signal decoder
% https://dl.acm.org/doi/10.1145/3546869
% https://github.com/jkadbear/LoRaPHY/blob/master/LoRaPHY.m
% https://github.com/tapparelj/gr-lora_sdr/blob/master/lib
%
% rx_sdr -g12 -f 910300000 -s 250000 -F CF32 /tmp/lora.raw
% phy = LoraPhy(7, 125e3, '/tmp/lora.raw', 1024e3)
% [x, netid1, netid2] = phy.detect_preamble();
% [sfd, hdr] = phy.detect_sfd(x);
% [payload_len, cr, crc, is_valid] = phy.decode_header(hdr);
% phy.plot_symbols(1024, 13)
%

classdef LoraPhy < handle & matlab.mixin.Copyable
    properties (SetAccess = immutable)
        sf;            % spreading factor (7,8,9,10,11,12)
        sr;            % sample rate: 2x
        bw;            % bandwidth 125 kHz
        fs;            % sample freq: sr * bw = 250 kHz
        cr;            % coding rate: (1:4/5 2:4/6 3:4/7 4:4/8)
        use_ldro;      % low data rate optimization
        use_crc;       % calculate cyclic redundancy check
        use_hamming;   % hamming error detection and correction
        has_header;    % data has header
        preamble_len;  % preamble length
    end

    properties (Access = private)
        sig;           % lora signal
        sps;           % samples per symbol
        ft_ratio;      % fft bins per sps
        ft_sym_fct     % fft symbol factor
        ft_bins;       % number of fft bins
        whitening_seq; % lfsr: x^8+x^6+x^5+x^4+1
        upchirp;       % downchirp removal
        downchirp;     % upchirp removal
    end

    methods
        function this = LoraPhy(sf, bw, filename, file_fs, swap_iq)
            % LoraPhy
            if(nargin < 5)
                swap_iq = false;
            end

            this.sf = sf;                   % spreading factor (7,8,9,10,11,12)
            this.sr = 2;                    % sample rate: 2x
            this.bw = bw;                   % bandwidth 125 kHz
            this.fs = this.sr * bw;         % sample freq: sr * bw = 250 kHz
            this.cr = 1;                    % coding rate: (1:4/5 2:4/6 3:4/7 4:4/8)
            this.use_ldro = 0;              % low data rate optimization (0/1)
            this.use_crc = 1;               % cyclic redundancy check (0/1)
            this.use_hamming = 1;           % hamming error detection and correction (0/1)
            this.has_header = 1;            % decode header (0/1)
            this.preamble_len = 8;          % preamble length

            this.sps = this.sr * 2^sf;      % samples per symbol
            this.ft_ratio = 2;              % fft bins per sps
            this.ft_sym_fct = this.sr * this.ft_ratio; % fft symbol factor
            this.ft_bins = this.ft_ratio * this.sps;   % fft bins per symbol

            % https://github.com/tapparelj/gr-lora_sdr/blob/master/lib/tables.h
            % linear-feedback shift register
            % x^8+x^6+x^5+x^4+1
            % reg = 0xff;
            % for ii = 1:255
            %     fprintf('0x%02x, ', reg);
            %     if(~mod(ii, 16))
            %         fprintf('\n');
            %     end
            %     reg = bitxor(bitshift(reg,1), bitxor(bitget(reg,4), bitxor(bitget(reg,5), bitxor(bitget(reg,6), bitget(reg,8)))));
            % end
            this.whitening_seq = uint8([ ...
                0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0x0b, 0x17, 0x2f, 0x5e, 0xbc, 0x78, 0xf1, 0xe3, ...
                0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47, ...
                0x8e, 0x1c, 0x38, 0x71, 0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0, ...
                0xc0, 0x81, 0x03, 0x06, 0x0c, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4d, 0x9b, ...
                0x37, 0x6e, 0xdc, 0xb9, 0x72, 0xe4, 0xc8, 0x90, 0x20, 0x41, 0x82, 0x05, 0x0a, 0x15, 0x2b, 0x56, ...
                0xad, 0x5b, 0xb6, 0x6d, 0xda, 0xb5, 0x6b, 0xd6, 0xac, 0x59, 0xb2, 0x65, 0xcb, 0x96, 0x2c, 0x58, ...
                0xb0, 0x61, 0xc3, 0x87, 0x0f, 0x1f, 0x3e, 0x7d, 0xfb, 0xf6, 0xed, 0xdb, 0xb7, 0x6f, 0xde, 0xbd, ...
                0x7a, 0xf5, 0xeb, 0xd7, 0xae, 0x5d, 0xba, 0x74, 0xe8, 0xd1, 0xa2, 0x44, 0x88, 0x10, 0x21, 0x43, ...
                0x86, 0x0d, 0x1b, 0x36, 0x6c, 0xd8, 0xb1, 0x63, 0xc7, 0x8f, 0x1e, 0x3c, 0x79, 0xf3, 0xe7, 0xce, ...
                0x9c, 0x39, 0x73, 0xe6, 0xcc, 0x98, 0x31, 0x62, 0xc5, 0x8b, 0x16, 0x2d, 0x5a, 0xb4, 0x69, 0xd2, ...
                0xa4, 0x48, 0x91, 0x22, 0x45, 0x8a, 0x14, 0x29, 0x52, 0xa5, 0x4a, 0x95, 0x2a, 0x54, 0xa9, 0x53, ...
                0xa7, 0x4e, 0x9d, 0x3b, 0x77, 0xee, 0xdd, 0xbb, 0x76, 0xec, 0xd9, 0xb3, 0x67, 0xcf, 0x9e, 0x3d, ...
                0x7b, 0xf7, 0xef, 0xdf, 0xbf, 0x7e, 0xfd, 0xfa, 0xf4, 0xe9, 0xd3, 0xa6, 0x4c, 0x99, 0x33, 0x66, ...
                0xcd, 0x9a, 0x35, 0x6a, 0xd4, 0xa8, 0x51, 0xa3, 0x46, 0x8c, 0x18, 0x30, 0x60, 0xc1, 0x83, 0x07, ...
                0x0e, 0x1d, 0x3a, 0x75, 0xea, 0xd5, 0xaa, 0x55, 0xab, 0x57, 0xaf, 0x5f, 0xbe, 0x7c, 0xf9, 0xf2, ...
                0xe5, 0xca, 0x94, 0x28, 0x50, 0xa1, 0x42, 0x84, 0x09, 0x13, 0x27, 0x4f, 0x9f, 0x3f, 0x7f, ]');

            this.sig = LoraPhy.read_cu8(filename);
            this.sig = lowpass(this.sig, this.bw/2, file_fs);
            if(this.fs ~= file_fs)
                this.sig = resample(this.sig, this.fs, file_fs);  % resample signal @ 2x bandwidth
            end
            if(swap_iq)
                % swap I and Q channels
                this.sig = imag(this.sig) + 1i * real(this.sig);
            end

            t = (0:this.sps-1) / this.fs;
            this.upchirp = chirp(t, -this.bw/2, t(end), this.bw/2, 'linear', 0, 'complex').';
            %close all; pspectrum(this.upchirp, this.fs, 'spectrogram', Reassign=true);
            this.downchirp = chirp(t, this.bw/2, t(end), -this.bw/2, 'linear', 0, 'complex').';
        end

        function [x, netid1, netid2] = detect_preamble(this, pos, preamble_len, invert)
            % DETECT_PREAMBLE(pos, [invert], [preamble_len])
            % input:
            %   pos:           offset from the beginning of signal
            %   invert:        detect downchirps rather than upchirps
            %   preamble_len:  preamble length
            % output:
            %    x:  signal position of fft bin0 for the last detected preamble chirp
            x = 0;
            netid1 = 0;
            netid2 = 0;

            if((nargin < 2) || (pos < 1))
                pos = 1;
            end
            if(nargin < 3)
                preamble_len = this.preamble_len;
            end
            if(nargin < 4)
                invert = false;
            end

            det_count = 1; % current sample == 1
            fbin_last = 0;
            while(pos <= (length(this.sig) - this.sps))
                [fbin, mval] = this.dechirp(pos, invert);
                %fprintf('%d) detect_preamble (%4d) - fbin: %3d  mval: %3d\n', det_count, pos, fbin, mval);

                dfbin = abs(fbin - fbin_last);
                if((mval > 15) && ((dfbin <= this.ft_ratio) || (dfbin >= (this.ft_bins-this.ft_ratio))))
                    det_count = det_count + 1;
                    pos_adj = mod(round((fbin-1)/this.ft_ratio), this.sps);
                    if(pos_adj > 16)
                        % move fbin to 1
                        pos = pos - pos_adj;
                        fbin = 1;
                    end

                    if(det_count >= preamble_len)
                        % preamble detected, adjust fft bin to begin at 1
                        pos_adj = mod(round((fbin-1)/this.ft_ratio), this.sps);
                        if(invert)
                            x = pos + fbin - 1;  % inverted chirp, non-inverted IQ
                        else
                            x = pos - pos_adj;  % non-inverted chirp, inverted IQ
                        end

                        % consume any extra preamble chirps while looking
                        % for network id 1
                        while((dfbin <= this.ft_ratio) || (dfbin >= (this.ft_bins-this.ft_ratio)))
                            x = x + this.sps;
                            fbin1 = this.dechirp(x, invert);
                            dfbin = abs(fbin1 - fbin_last);
                        end

                        x = x + this.sps;
                        fbin2 = this.dechirp(x, invert);

                        if(invert)
                            fbin1 = this.sps - fbin1 + 1;  % inverted chirp, non-inverted IQ
                            fbin2 = this.sps - fbin2 + 1;
                        else
                            fbin1 = fbin1 - 1;             % non-inverted chirp, inverted IQ
                            fbin2 = fbin2 - 1;
                        end

                        netid1 = round(fbin1 / this.ft_sym_fct);
                        netid2 = round(fbin2 / this.ft_sym_fct);

                        %fprintf('  *** preamble detected ***  fbin:%4d  x:%d  netid1:%d  netid2:%d\n', fbin, x, netid1, netid2);
                        return;
                    end
                else
                    det_count = 1;
                end

                fbin_last = fbin;
                pos = pos + this.sps;
            end
        end

        % look for the start of frame delimiter
        function [sfd, hdr] = detect_sfd(this, pos, invert)
            sfd = 0;
            hdr = 0;

            if(pos < 1)
                return;
            end

            if(nargin < 3)
                invert = false;
            end

            while(pos <= (length(this.sig) - this.sps))
                [~, mval_up] = this.dechirp(pos, invert);
                [~, mval_dn] = this.dechirp(pos, ~invert);
                if mval_dn > mval_up
                    % downchirp detected
                    % look for second downchirp at 1.25 sps
                    [~, mval_up] = this.dechirp(pos + 1.25 * this.sps, invert);
                    [~, mval_dn] = this.dechirp(pos + 1.25 * this.sps, ~invert);
                    if mval_dn > mval_up
                        % second downchirp detected
                        sfd = pos;
                        hdr = sfd + 2.25 * this.sps;  % the sfd is 2.25 symbols long
                        return;
                    end
                end

                pos = pos + this.sps;
            end
        end

        function [payload_len, cr, crc, is_valid] = decode_header(this, pos, invert)
            if(pos < 1)
                return;
            end

            if(nargin < 3)
                invert = false;
            end

            symbols = zeros(8, 1);
            for ii = 1:8
                fbin = this.dechirp(pos, invert);

                if(invert)
                    fbin = this.sps - fbin + 1;  % inverted chirp, non-inverted IQ
                else
                    fbin = fbin - 1;             % non-inverted chirp, inverted IQ
                end

                symbols(ii) = round(fbin / this.ft_sym_fct);
                %fprintf('%d) pos:%4d  fbin:%3d  -->  sym:%3d\n', ii, pos, fbin, symbols(ii));

                pos = pos + this.sps;
            end

            % gray decoding
            symbols_g = this.gray_decode(symbols, true);

            % deinterleave
            codewords = this.diag_deinterleave(symbols_g, 4, this.sf-2);

            % hamming decode
            header = this.hamming_decode(codewords, 4);

            % parse header
            payload_len = bitor(bitshift(header(1), 4), header(2));
            crc = bitand(header(3), 1);
            cr = bitshift(header(3), -1);

            % validate the header checksum
            header_checksum = [bitand(header(4), 1); int2bit(header(5), 4)];
            header_checksum_calc = this.calc_header_csum(header);
            is_valid = all(header_checksum == header_checksum_calc);
        end

        function symbols = decode_payload(this, pos_hdr, payload_len, invert)
            if(nargin < 4)
                invert = false;
            end

            symcnt = this.calc_payload_symbol_count(payload_len);

            pos = pos_hdr + 8*this.sps;
            symbols = zeros(symcnt,1,'uint16');
            for ii = 1:symcnt
                if(pos > (length(this.sig) - this.sps))
                    fprintf('error: unexpected end of data reached - pos:%d  sym_num:%d\n', pos, ii);
                    return;
                end

                fbin = this.dechirp(pos, invert);

                if(invert)
                    fbin = this.sps - fbin + 1;  % inverted chirp, non-inverted IQ
                else
                    fbin = fbin - 1;             % non-inverted chirp, inverted IQ
                end

                symbols(ii) = round(fbin / this.ft_sym_fct);
                %fprintf('%d) pos:%4d  fbin:%3d  -->  sym:%3d\n', ii, pos, fbin, symbols(ii));

                pos = pos + this.sps;
            end
        end

        function symcnt = calc_payload_symbol_count(this, payload_len)
            % calc_payload_symbols  calculate number of symbols
            %
            % input:
            %     payload_len: payload length
            % output:
            %     symcnt: number of symbols required to encode payload
            %     length

            % see SX1276_Datasheet.pdf p.31
            symcnt = max(ceil((2*double(payload_len)-this.sf +7+ 4*this.use_crc-5*(1-this.has_header)) / (this.sf-2*this.use_ldro)) * (this.cr+4), 0);
        end

        function hcsum = calc_header_csum(this, header)
            hdata = int2bit(header(1:3)', 4, false)';

            hcsum = zeros(5,1,'uint8');
            hcsum(1) = this.freduce(@xor, [hdata(1,4), hdata(1,3), hdata(1,2), hdata(1,1)]);
            hcsum(2) = this.freduce(@xor, [hdata(1,4), hdata(2,4), hdata(2,3), hdata(2,2), hdata(3,1)]);
            hcsum(3) = this.freduce(@xor, [hdata(1,3), hdata(2,4), hdata(2,1), hdata(3,4), hdata(3,2)]);
            hcsum(4) = this.freduce(@xor, [hdata(1,2), hdata(2,3), hdata(2,1), hdata(3,3), hdata(3,2), hdata(3,1)]);
            hcsum(5) = this.freduce(@xor, [hdata(1,1), hdata(2,2), hdata(3,4), hdata(3,3), hdata(3,2), hdata(3,1)]);
        end

        function symbols_g = gray_decode(this, symbols, use_ldro)
            if(use_ldro)
                sym_b2 = bitshift(uint16(symbols), -2);
            else
                sym_b2 = uint16(mod(int16(symbols)-1, 2^this.sf));
            end
            sym_b3 = bitshift(sym_b2, -1);
            symbols_g = bitxor(sym_b2, sym_b3);
        end

        function codewords = diag_deinterleave(this, symbols_g, cr, sf)
            % DIAG_DEINTERLEAVE(symbols_g, bits)
            %   perform circular left shift by n bits
            %                                           154 0 163 92 0
            %   |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
            %    0 |1  0  1  0 : 0| 1  0  1  0  <<1  ->  1  0  1  0  0
            %    1  0 |0  0  0 : 1  0| 0  0  0  <<2  ->  0  0  0  1  0
            %    0  1  0 |1  0 : 0  1  0| 1  0  <<3  ->  1  0  0  1  0
            %    0  0  1  0 |1 : 0  0  1  0| 1  <<4  ->  1  0  0  1  0
            %   |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
            %    0 |0  0  0  1 : 0| 0  0  0  1  <<1  ->  0  0  0  1  0
            %    0  0 |1  0  1 ; 0  0| 1  0  1  <<2  ->  1  0  1  0  0
            if(nargin < 3)
                cr = this.cr;
            end
            if(nargin < 4)
                sf = this.sf;
            end

            cr_bits = cr + 4;
            cw_len = (length(symbols_g) / cr_bits) * sf;
            codewords = zeros(cw_len,1,'uint8');

            cwi = 1;
            for offs = 1:cr_bits:length(symbols_g)
                temp = zeros(cr_bits,sf,'uint8');
                for ii = (1:cr_bits)
                    temp(ii,:) = circshift(int2bit(symbols_g(offs+cr_bits-ii), sf, false)', cr_bits-ii);
                end
                codewords(cwi:cwi+sf-1,1) = bit2int(temp, cr_bits)';
                cwi = cwi + sf;
            end
        end

        function data = hamming_decode(this, codewords, cr)
            % HAMMING_DECODE(codewords, cr)
            %
            %           parity    data
            %  cr 4/8  p3p2p1p0 d3d2d1d0
            %  cr 4/7    p2p1p0 d3d2d1d0
            %  cr 4/6      p1p0 d3d2d1d0
            %  cr 4/5        p4 d3d2d1d0
            %
            %    p0 = d0 ^ d1 ^ d2
            %    p1 = d1 ^ d2 ^ d3
            %    p2 = d0 ^ d1 ^ d3
            %    p3 = d0 ^ d2 ^ d3
            %    p4 = d0 ^ d1 ^ d2 ^ d3
            if(nargin < 3)
                cr = this.cr;
            end
            cr_bits = cr + 4;

            data = uint8(bitand(uint8(codewords), 0x0f));
            if(~this.use_hamming)
                return;
            end

            % TODO: support CR 4/7, 4/6, and 4/5 (not needed to decode header)
            % p3p2p1p0 -> p1p2p5p3
            if(cr_bits < 8)
                return;
            end

            parity = bitand(bitshift(codewords, -4), 0x0f);

            % calculate parity
            d0 = bitand(data, 1);
            d1 = bitand(bitshift(data, -1), 1);
            d2 = bitand(bitshift(data, -2), 1);
            d3 = bitand(bitshift(data, -3), 1);

            p0 = bitxor(d0, bitxor(d1, d2));
            p1 = bitxor(d1, bitxor(d2, d3));
            p2 = bitxor(d0, bitxor(d1, d3));
            p3 = bitxor(d0, bitxor(d2, d3));
            %p4 = bitxor(bitxor(d0, d1), bitxor(d2, d3));  % CR 4/5

            pcalc = bitor(bitor(bitshift(p3, 3), bitshift(p2, 2)), bitor(bitshift(p1, 1), p0));

            % check / repair data
            for ii = 1:length(data)
                perr = bitxor(parity(ii), pcalc(ii));
                if(perr)
                    % repair data
                    perr = bitxor(perr, 0x0f);
                    % bit 1->8, bit 8->2, bit 2->1
                    be = bitor( ...
                        bitor(bitshift(bitand(perr,0x01),  3),          bitand(perr,0x04)), ...
                        bitor(bitshift(bitand(perr,0x08), -2), bitshift(bitand(perr,0x02), -1)) );
                    if(ismember(be, [1, 2, 4, 8]))
                        data(ii) = bitxor(data(ii), be);
                        fprintf('parity data correction - data elem:%d  bit:%d  repaired data:%d\n', ii, be, data(ii));
                    else
                        fprintf(' !!!! unrecoverable parity error !!!! - data elem:%d  bit:%d\n', ii, be);
                    end
                end
            end
        end

        function bytes_w = dewhiten(self, bytes)
            % DEWHITEN(bytes)  Data Dewhitening
            %
            % input:
            %     bytes: Bytes after deinterleaving
            % output:
            %     bytes_w: Bytes after dewhitening
            len = length(bytes);
            bytes_w = bitxor(uint8(bytes(1:len)), self.whitening_seq(1:len));
        end

        function [fbin, mval, ft_pow] = dechirp(this, pos, invert)
            % DECHIRP(pos, [invert])
            % apply dechirping on the symbol at pos
            % input:
            %      pos:  starting signal position
            %   invert:  apply upchirp rather than downchirp
            % output:
            %     fbin:  fft bin with the greatest power
            %     mval:  maximum fft power value
            %   ft_pow:  vector of power levels by bin
            if(nargin < 3)
                invert = false;
            end

            if(invert)
                chp = this.upchirp;
            else
                chp = this.downchirp;
            end

            ft = fft(this.sig(pos:pos+this.sps-1) .* chp, 2*this.ft_bins);
            ft_pow = abs(ft(1:this.ft_bins)) + abs(ft(this.ft_bins+1:2*this.ft_bins));

            [mval, fbin] = max(ft_pow);
            mval = uint16(round(mval));
        end

        function plot_symbols(this, pos, sym_count, use_legend)
            % PLOT_SYMBOLS(this, pos, [sym_count])
            %   pos:        offset from the beginning of signal
            %   sym_count:  the number of symbols to plot (default 9)
            %
            %   create upchirped / downcirped plots of the signal starting at
            %    'pos' for 'sym_count' symbols
            %   note: the max fft power value is indicated by a red diamond
            if(nargin < 3)
                sym_count = 9;
            end
            if(nargin < 4)
                use_legend = false;
            end

            close(findobj('type', 'figure', 'number', 1));
            f = figure(1);
            scn = get(0, 'screensize');
            f.Position = [20,400,scn(1,3)-120,680];
            tiledlayout(2, sym_count, 'TileSpacing', 'compact', 'Padding', 'compact');

            axes = zeros(2*sym_count,1);
            type = ["up-chirped ", "down-chirped "];
            for cdir = 1:2
                tpos = pos;
                for jj = 1:sym_count
                    [fbin, mval, ft_pow] = this.dechirp(tpos, (cdir==1));

                    ax = nexttile;
                    plot(ft_pow);
                    %hist(ft_pow, length(ft_pow));
                    title(type(cdir) + jj);
                    ax.XLim = [0, this.ft_bins];
                    ax.XTick = linspace(0, this.ft_bins, 5);
                    axes((cdir-1)*sym_count + jj) = ax;

                    hold on
                    mval = uint16(mval);
                    plot(fbin, mval, 'rd');
                    if(use_legend)
                        legend('', sprintf('%d,%d', fbin, mval));
                    else
                        if(fbin > this.sps/2)
                            txtx = fbin - 120;
                        else
                            txtx = fbin + 20;
                        end
                        text(txtx, mval, sprintf('%d,%d', fbin, mval));
                    end
                    %fprintf('%4d) plot(%d,%2d)  fbin:%3d  sym:%3d\n', tpos, cdir, jj, fbin, uint16((fbin-1)/2));

                    tpos = tpos + this.sps;
                end
            end

            linkaxes(axes, 'xy');
        end

        function plot_symbol_pspec(this, pos, count)
            if(nargin < 3)
                count = 1;
            end

            pspectrum(this.sig(pos:pos+this.sps*count), this.fs, 'spectrogram', Reassign=true);
        end
    end

    methods(Static)
        function v = freduce(f, a)
            % FREDUCE(f, a)
            %   reduce the array 'a' using the function 'f'
            v = a(1);
            for ii = 2:length(a)
                v = f(v, a(ii));
            end
        end

        function v = read_cu8(filename, count)
            % READ(filename, [count])
            %   open filename and return the contents as a column vector,
            %   treating them as 32 bit complex numbers
            narginchk(1, 2);

            if(nargin < 2)
                count = Inf;
            end

            f = fopen(filename, 'rb');
            if(f < 0)
                v = 0;
            else
                t = fread(f, [2, count], 'uint8');
                fclose(f);
                t = (single(t) - 127.0) / 128.0;
                v = t(1,:) + t(2,:)*1i;
                [r, c] = size(v);
                v = reshape(v, c, r);
            end
        end

        function v = read_cf32(filename, count)
            % READ(filename, [count])
            %   open filename and return the contents as a column vector,
            %   treating them as 32 bit complex numbers
            %   https://github.com/gnuradio/gnuradio/blob/master/gr-utils/octave/read_complex_binary.m
            narginchk(1, 2);

            if(nargin < 2)
                count = Inf;
            end

            f = fopen(filename, 'rb');
            if(f < 0)
                v = 0;
            else
                t = fread(f, [2, count], 'float32');
                fclose(f);
                v = t(1,:) + t(2,:)*1i;
                [r, c] = size(v);
                v = reshape(v, c, r);
            end
        end

        function v = write_cf32(filename, data)
            % WRITE(filename, data)
            %   write IQ data to 'filename'
            %   data format is interleaved 32-bit IQ float values: IQIQIQ....
            %   note: this is compatible with READ()
            %   https://github.com/gnuradio/gnuradio/blob/master/gr-utils/octave/write_complex_binary.m
            narginchk(2, 2);

            f = fopen(filename, 'wb');
            if(f < 0)
                v = 0;
            else
                re = real(data);
                im = imag(data);
                re = re(:)';
                im = im(:)';
                y = [re;im];
                y = y(:);
                v = fwrite(f, y, 'float');
                fclose(f);
            end
        end
    end
end

