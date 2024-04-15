
% lora signal decoder
% https://dl.acm.org/doi/10.1145/3546869
% https://github.com/jkadbear/LoRaPHY/blob/master/LoRaPHY.m
% https://github.com/tapparelj/gr-lora_sdr/blob/master/lib

% phy = LoraPhy(7, 125e3, 'lora_923.3_sample/lora.raw', 1024e3)
% [x, netid1, netid2] = phy.detect_preamble();
% [sfd, hdr] = phy.detect_sfd(x);
% [payload_len, cr, crc, is_valid] = phy.decode_header(hdr);
% phy.plot_symbols(1024, 13)

classdef LoraPhy < handle & matlab.mixin.Copyable
    properties (SetAccess = immutable)
        sf;            % spreading factor (7,8,9,10,11,12)
        bw;            % bandwidth 125 kHz
        fs;            % sample rate: 2x bw = 250 kHz
        cr;            % coding rate: (1:4/5 2:4/6 3:4/7 4:4/8)
        preamble_len;  % preamble length
        ft_det_bins;   % num fft detect bins
    end

    properties  (Access = private)
        sig;        % lora signal
        sps;        % samples per symbol
        ft_bins;    % number of fft bins
        ft_len;     % fft size
        cfo;        % carrier frequency offset
        upchirp;    % downchirp removal
        downchirp;  % upchirp removal
    end

    methods
        function this = LoraPhy(sf, bw, filename, file_fs, swap_iq)
            % LoraPhy
            if(nargin < 5)
                swap_iq = false;
            end

            this.sf = sf;             % spreading factor (7,8,9,10,11,12)
            this.bw = bw;             % bandwidth 125 kHz
            this.fs = bw * 2;         % sample rate: 2x bw = 250 kHz
            this.cr = 1;              % coding rate: (1:4/5 2:4/6 3:4/7 4:4/8)
            this.preamble_len = 8;    % preamble length
            this.ft_det_bins = 4;     % num fft detect bins

            this.sps = 2 * 2^sf;      % samples per symbol
            this.ft_bins = this.sps;  % number of fft bins
            this.ft_len = 2 * this.ft_bins;  % fft size
            this.cfo = 0;             % carrier frequency offset

            this.sig = LoraPhy.read(filename);
            this.sig = lowpass(this.sig, bw/2, file_fs);
            this.sig = resample(this.sig, this.fs, file_fs);  % resample signal @ 2x bandwidth
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
                fbin = this.dechirp(pos, invert);
                %fprintf('%4d) detect_preamble - fbin:%d\n', pos, fbin);

                if(abs(fbin - fbin_last) <= this.ft_det_bins)
                    det_count = det_count + 1;
                    if(det_count >= preamble_len)
                        % preamble detected, adjust fft bin to zero
                        if(invert)
                            x = pos + fbin - 1;             % inverted chirp, non-inverted IQ
                        else
                            x = pos + this.sps - fbin + 1;  % non-inverted chirp, inverted IQ
                        end

                        % read network id
                        fbin1 = this.dechirp(x+this.sps, invert);
                        fbin2 = this.dechirp(x+2*this.sps, invert);
                        if(invert)
                            fbin1 = this.sps - fbin1 + 1;  % inverted chirp, non-inverted IQ
                            fbin2 = this.sps - fbin2 + 1;
                        else
                            fbin1 = fbin1 - 1;             % non-inverted chirp, inverted IQ
                            fbin2 = fbin2 - 1;
                        end

                        netid1 = fbin1 / 2;
                        netid2 = fbin2 / 2;

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

                symbols(ii) = fbin / 2;
                %fprintf('%d) pos:%4d  fbin:%3d  -->  sym:%3d\n', ii, pos, fbin, symbols(ii));

                pos = pos + this.sps;
            end

            % gray decoding
            symbols_g = this.gray_decode(symbols);

            % deinterleave
            codewords = this.diag_deinterleave(symbols_g, this.sf-2);

            % hamming decode
            header = this.hamming_decode(codewords, 8);

            % parse header
            payload_len = bitor(bitshift(header(1), 4), header(2));
            crc = bitand(header(3), 1);
            cr = bitshift(header(3), -1);

            % validate the header checksum
            header_checksum = [bitand(header(4), 1); int2bit(header(5), 4)];
            header_checksum_calc = this.calc_header_csum(header);
            is_valid = all(header_checksum == header_checksum_calc);
        end

        function hcsum = calc_header_csum(this, header)
            hdata = int2bit(header(1:3)', 4, false)';

            hcsum = zeros(5,1,'uint8');
            hcsum(5) = this.freduce(@xor, [hdata(1,4), hdata(1,3), hdata(1,2), hdata(1,1)]);
            hcsum(4) = this.freduce(@xor, [hdata(1,4), hdata(2,4), hdata(2,3), hdata(2,2), hdata(3,1)]);
            hcsum(3) = this.freduce(@xor, [hdata(1,3), hdata(2,4), hdata(2,1), hdata(3,4), hdata(3,2)]);
            hcsum(2) = this.freduce(@xor, [hdata(1,2), hdata(2,3), hdata(2,1), hdata(3,3), hdata(3,2), hdata(3,1)]);
            hcsum(1) = this.freduce(@xor, [hdata(1,1), hdata(2,2), hdata(3,4), hdata(3,3), hdata(3,2), hdata(3,1)]);
        end

        function symbols_g = gray_decode(this, symbols)
            sym_b2 = bitshift(uint16(symbols), -2);
            sym_b3 = bitshift(sym_b2, -1);
            symbols_g = bitxor(sym_b2, sym_b3);
        end

        function codewords = diag_deinterleave(this, symbols_g, bits)
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
            codewords = bit2int([
                circshift(int2bit(symbols_g(8), bits, false)', 7);
                circshift(int2bit(symbols_g(7), bits, false)', 6);
                circshift(int2bit(symbols_g(6), bits, false)', 5);
                circshift(int2bit(symbols_g(5), bits, false)', 4);
                circshift(int2bit(symbols_g(4), bits, false)', 3);
                circshift(int2bit(symbols_g(3), bits, false)', 2);
                circshift(int2bit(symbols_g(2), bits, false)', 1);
                circshift(int2bit(symbols_g(1), bits, false)', 0);
            ], 8)';
        end
        % function codewords = diag_deinterleave(this, symbols_g, bits)
        %     temp = zeros(8,1);
        %     for ii = (8:1)
        %         temp(ii,:) = circshift(int2bit(symbols_g(ii), bits, false)', ii-1);
        %     end
        %     codewords = bit2int(temp,8);
        % end

        % TODO: support CR 4/7, 4/6, and 4/5 (not needed to decode header)
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
            data = bitand(codewords, 0x0f);
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
                        fprintf("parity data correction - data elem:%d  bit:%d  repaired data:%d\n", ii, be, data(ii));
                    else
                        fprintf(" !!!! unrecoverable parity error !!!! - data elem:%d  bit:%d\n", ii, be);
                    end
                end
            end
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

            ft = fft(this.sig(pos:pos+this.sps-1) .* chp, this.ft_len);
            ft_pow = abs(ft(1:this.ft_bins)) + abs(ft(this.ft_len - this.ft_bins+1:this.ft_len));

            [mval, fbin] = max(ft_pow);
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

            axes = zeros(2 * sym_count, 1);
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
                    mval = int8(mval);
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
                    %fprintf('%4d) plot(%d,%2d)  fbin:%3d  sym:%3d\n', tpos, cdir, jj, fbin, int8((fbin-1)/2));

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

        function v = read(filename, count)
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
                t = fread(f, [2, count], 'float');
                fclose(f);
                v = t(1,:) + t(2,:)*1i;
                [r, c] = size(v);
                v = reshape(v, c, r);
            end
        end

        function v = write(filename, data)
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
