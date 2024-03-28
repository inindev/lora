
% a copy of: https://github.com/jkadbear/LoRaPHY/blob/master/LoRaPHY.m
% rewritten to simplify and make it easier for me to read
%

% phy = LoraPhy(7, 125e3, 'lora_923.3_sample/lora.raw', 1024e3, true)
% phy.detect(1)
% for jj = (1:256); phy.detect(jj); end
% phy.plot_symbols(1024, 13)

classdef LoraPhy < handle & matlab.mixin.Copyable
    properties (SetAccess = immutable)
        sf;            % spreading factor (7,8,9,10,11,12)
        bw;            % bandwidth 125 kHz
        fs;            % sample rate: 2x bw = 250 kHz
        cr;            % code rate: (1:4/5 2:4/6 3:4/7 4:4/8)
        preamble_len;  % preamble length
        ft_det_bins;   % num fft detect bins
    end

    properties  (Access = private)
        sig;        % lora signal
        sps;        % samples per symbol
        ft_bins;    % number of fft bins
        ft_len;     % fft size
        cfo;        % carrier frequency offset
        upchirp;    % used for de-chirping
        downchirp;  % used for de-chirping
        chirps;     % used for de-chirping
    end

    methods
        function this = LoraPhy(sf, bw, filename, file_fs, invert)
            % LoraPhy
            if(nargin < 5)
                invert = false;
            end

            this.sf = sf;             % spreading factor (7,8,9,10,11,12)
            this.bw = bw;             % bandwidth 125 kHz
            this.fs = bw * 2;         % sample rate: 2x bw = 250 kHz
            this.cr = 1;              % code rate: (1:4/5 2:4/6 3:4/7 4:4/8)
            this.preamble_len = 8;    % preamble length
            this.ft_det_bins = 4;     % num fft detect bins

            this.sps = 2 * 2^sf;      % samples per symbol
            this.ft_bins = this.sps;  % number of fft bins
            this.ft_len = 2 * this.ft_bins;  % fft size
            this.cfo = 0;             % carrier frequency offset

            this.sig = LoraPhy.read(filename);
            this.sig = lowpass(this.sig, bw/2, file_fs);
            this.sig = resample(this.sig, this.fs, file_fs);  % resample signal @ 2x bandwidth
            if(invert)
                % swap I and Q channels
                this.sig = imag(this.sig) + 1j * real(this.sig);
            end

            t = (0:this.sps-1) / this.fs;
            this.upchirp = chirp(t, -this.bw/2, t(end), this.bw/2, 'linear', 0, 'complex').';
            %close all; pspectrum(this.upchirp, this.fs, "spectrogram", Reassign=true);
            this.downchirp = chirp(t, this.bw/2, t(end), -this.bw/2, 'linear', 0, 'complex').';
        end

        function x = detect(this, pos)
            % DETECT(pos, cdir)
            %   pos:    offset from the beginning of signal
            x = 0;
            if(pos < 1)
                return;
            end

            start = pos;
            det_count = 0;
            last = 0;
            while(pos <= (length(this.sig) - this.sps))
                fbin = dechirp(this, pos);

                if(abs(fbin - last) <= this.ft_det_bins)
                    det_count = det_count + 1;
                    if(det_count > this.preamble_len-2)
                        % preamble detected, adjust fft bin to zero
                        x = pos - fbin + 1;
                        fprintf("start:%3d  fbin:%4d  x:%d  *** preamble detected ***\n", start, fbin, x);
                        return;
                    end
                else
                    det_count = 0;
                end
                last = fbin;
                pos = pos + this.sps;
            end
        end

        function [fbin, mval, ft_pow] = dechirp(this, pos, is_up)
            % DECHIRP:  apply dechirping on the symbol at pos
            % input:
            %      pos:  symbol index
            %    is_up:  apply upchirp rather than downchirp
            % output:
            %     fbin:  fft bin with the greatest power
            %     mval:  maximum fft power value
            %   ft_pow:  vector of power levels by bin
            if(nargin < 3)
                is_up = false;
            end

            if(is_up)
                chp = this.upchirp;
            else
                chp = this.downchirp;
            end

            ft = fft(this.sig(pos:pos+this.sps-1) .* chp, this.ft_len);
            ft_pow = abs(ft(1:this.ft_bins)) + abs(ft(this.ft_len - this.ft_bins+1:this.ft_len));
            [mval, fbin] = max(ft_pow);
        end

        function plot_symbols(this, pos, sym_count)
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

            type = ["upchirped", "dnchirped"];

            f = figure(1);
            scn = get(0,'screensize');
            f.Position = [1,439,scn(1,3)-120,680];

            for jj = 1:sym_count
                for cdir = 1:2
                    % 2 x n plot
                    [fbin, mval, ft_pow] = this.dechirp(pos, (cdir==1));
                    subplot(2, sym_count, ((cdir-1) * sym_count + jj));
                    plot(ft_pow);
                    %hist(ft_pow, length(ft_pow));
                    title(type(cdir) + " " + jj);

                    hold on
                    plot(fbin, mval, 'rd');
                    %fprintf("plot(%d,%d) fbin = %d\n", cdir, jj, fbin);
                end

                pos = pos + this.sps;
            end
        end
    end


    methods(Static)
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

