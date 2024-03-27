%
% a copy of: https://github.com/jkadbear/LoRaPHY/blob/master/LoRaPHY.m
% rewritten to simplify and make it easier for me to read
%

% phy = LoraPhy(7, 125e3, 1024e3, 'lora_923.3_sample/lora.raw');
% phy.detect(1, 1);
% for jj = (1:256); phy.detect(jj, 1); end
% phy.plot_symbols(1024, 13);

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
        function this = LoraPhy(sf, bw, fs, filename, invert)
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
            this.sig = lowpass(this.sig, bw/2, fs);
            this.sig = resample(this.sig, this.fs, fs);  % resample signal @ 2x bandwidth
            if(invert)
                % swap I and Q channels
                this.sig = imag(this.sig) + 1j * real(this.sig);
            end

            t = (0:this.sps-1) / this.fs;
            this.upchirp = chirp(t, -this.bw/2, t(end), this.bw/2, 'linear', 0, 'complex').';
            this.downchirp = chirp(t, this.bw/2, t(end), -this.bw/2, 'linear', 0, 'complex').';
            this.chirps = [this.upchirp, this.downchirp];
            %close all; pspectrum(this.chirps(:,1), this.fs, "spectrogram", Reassign=true);
        end

        %
        % dir: 0 = upchirp, 1 = downchirp
        function x = detect(this, pos, dir)
            x = 0;
            if(pos < 1)
                return;
            end

            start = pos;
            det_count = 0;
            last = 0;
            while(pos <= (length(this.sig) - this.sps))
                ft = fft(this.sig(pos:pos+this.sps-1) .* this.chirps(:,dir+1), this.ft_len);
                ft_pow = abs(ft(1:this.ft_bins)) + abs(ft(this.ft_len - this.ft_bins+1:this.ft_len));
                [~, fbin] = max(ft_pow);

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

        %
        % usage: plot_symbols(this, pos, [sym_count])
        %
        %  create upchirped / downcirped plots of the signal starting at
        %    'pos' for 'sym_count' symbols
        %  note: the max fft power value is indicated by a red diamond
        %
        function plot_symbols(this, pos, sym_count)
            if(nargin < 2)
                sym_count = 9;
            end

            type = ["upchirped", "dnchirped"];

            figure;

            for jj = 1:sym_count
                for dir = 0:1
                    % 2 x n plot
                    ft = fft(this.sig(pos:pos+this.sps-1) .* this.chirps(:,dir+1), this.ft_len);
                    ft_pow = abs(ft(1:this.ft_bins)) + abs(ft(this.ft_len - this.ft_bins+1:this.ft_len));

                    s = subplot(2, sym_count, dir*sym_count + jj);
                    plot(ft_pow);
                    %hist(ft_pow, length(ft_pow));
                    title(type(dir+1) + " " + jj);

                    hold on
                    [mval, fbin] = max(ft_pow);
                    plot(fbin, mval, 'rd');
                    %fprintf("plot(%d,%d) fbin = %d\n", dir+1, jj, fbin);
                end

                pos = pos + this.sps;
            end
        end

    end


    methods(Static)
        %
        % usage: read(filename, [count])
        %
        %  open filename and return the contents as a column vector,
        %  treating them as 32 bit complex numbers
        %
        % https://github.com/gnuradio/gnuradio/blob/master/gr-utils/octave/read_complex_binary.m
        function v = read(filename, count)
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

        %
        % usage: write(data, filename)
        %
        %  open filename and write data to it
        %  Format is interleaved float IQ e.g. each
        %  I,Q 32-bit float IQIQIQ....
        %  This is compatible with read()
        %
        % https://github.com/gnuradio/gnuradio/blob/master/gr-utils/octave/write_complex_binary.m
        function v = write(data, filename)
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

