%filename = 'motor_no_beads_1V_set.txt';
%filename = 'manual_drive_no_beads.txt';
%magnetic_spec = spectrogram_mag(filename, 3, 1, 'comb', true, true);
%plotspectrogram(magnetic_spec,50);

%%
%filename = 'poly_rot_M.txt'; %2 inside 1 outside (columns 6 and 3)
%filename = 'MGS1_air.txt'; %2 inside, 1 outide (columns 6 and 3)
%filename = '20s_field_mill_on.txt';%2 inside, 1 outside (columns 2 and 1)
%filename = "mgs1_drop_KM.txt";%2 inside, 1 outside (columns 2 and 1)
%filename = "MGS1C_rot.txt";
%filename = "MGS1C_FM_on.txt";
filename = "ljm_brass.csv";

%filename = 'testing.txt';%2 inside, 1 outside (columns 2 and 1)

magnetic_spec = spectrogram_mag2(filename, 2, 1, 'reduced', true, true);

figure
plotspectrogram(magnetic_spec,100);
%title('raw')

removed_spec = spectrogram_mag(filename, [2,1], 0, 'reduced', true, true);
removed_spec.t = removed_spec.t - 3.7; %mgs1
%removed_spec.t = removed_spec.t - 7.8; %mgs1c
figure(2)
plotspectrogram(removed_spec,100)

A=1;
function output = spectrogram_mag(filename, col, sc_flag, notch_mode, apply_highpass, apply_notch)
    % Defaults for optional arguments
    if nargin < 4 || isempty(notch_mode), notch_mode = 'comb'; end
    if nargin < 5 || isempty(apply_highpass), apply_highpass = true; end
    if nargin < 6 || isempty(apply_notch), apply_notch = true; end

    Fs = readmatrix(filename, "Range", "B2:B2");
    if isnan(Fs)
        Fs = readmatrix(filename, "Range", "B2:B2", 'LineEnding', 'Hz');
    end

    import = readmatrix(filename);
    if length(col) == 2
        data_1 = import(:, col(1)) - mean(import(:, col(1)));
        data_2 = import(:, col(2)) - mean(import(:, col(2)));
        data = (data_1 - data_2) * 1000;
    else
        data = import(:, col) * 1000;
        data = data - mean(data);
    end
    clear import;
    
    % High-pass filter
    if apply_highpass
        data = highpass_filter(data, Fs, 2);
    end

    % Notch filter
    if apply_notch
        switch lower(notch_mode)
            case 'butter'
                data = comb_notch_filter(data, Fs, 50, 5, 1);
            case 'comb'
                data = notch_50Hz(data, Fs);
            case 'reduced'
                Hd = notch_50_100(Fs, 1);
                data = Hd(data);
            otherwise
                warning('Unknown notch mode. Skipping notch filtering.');
        end
    end


    % Spectrogram settings
    M = Fs;
    L = Fs / 2;

    g = tukeywin(M,0.1);
    Ndft = 2^nextpow2(M);
    Ts = 1 / Fs;
    
    fft_single(data,Fs);

    [s, f, t] = spectrogram(data, g, L, Ndft, Fs, 'onesided', 'yaxis');
    s = s ./ M;
    s(2:end-1, :) = 2 * s(2:end-1, :);

    if sc_flag == 1
        load lemi_tf_dat.mat lemi_tf
        fd = (1 / (2 * pi)) * (2 / Ts) * atan((2 * pi * f * Ts) / 2);
        gain = interp1(lemi_tf(:, 1), lemi_tf(:, 2), fd, 'linear', 'extrap');
        phase = interp1(lemi_tf(:, 1), lemi_tf(:, 3), fd, 'linear', 'extrap');
        TF = (gain .* exp(1i * deg2rad(phase)))';
        sc = s ./ (ones(length(t), 1) * TF)';
    else
        sc = s .* 35;
    end

    output.s = sc;
    output.f = f;
    output.t = t;
end

function output = spectrogram_mag2(filename, col, sc_flag, notch_mode, apply_highpass, apply_notch)
    % Defaults for optional arguments
    if nargin < 4 || isempty(notch_mode), notch_mode = 'comb'; end
    if nargin < 5 || isempty(apply_highpass), apply_highpass = true; end
    if nargin < 6 || isempty(apply_notch), apply_notch = true; end
    
    fid = fopen(filename, 'r');
    Fs = fgetl(fid);
    fclose(fid);

    Fs = textscan(Fs, 'Sample rate: %f S/s');
    Fs = floor(Fs{1});

    import = readmatrix(filename);
    if length(col) == 2
        data_1 = import(:, col(1)) - mean(import(:, col(1)));
        data_2 = import(:, col(2)) - mean(import(:, col(2)));
        data = (data_1 - data_2) * 1000;
    else
        data = import(:, col) * 1000;
        data = data - mean(data);
    end
    clear import;
    
    % High-pass filter
    if apply_highpass
        data = highpass_filter(data, Fs, 2);
    end

    % Notch filter
    if apply_notch
        switch lower(notch_mode)
            case 'butter'
                data = comb_notch_filter(data, Fs, 50, 5, 1);
            case 'comb'
                data = notch_50Hz(data, Fs);
            case 'reduced'
                Hd = notch_50_100(Fs, 1);
                data = Hd(data);
            otherwise
                warning('Unknown notch mode. Skipping notch filtering.');
        end
    end


    % Spectrogram settings
    M = Fs;
    L = Fs / 2;

    L = 500;

    g = tukeywin(M,0.1);
    Ndft = 2^nextpow2(M);
    Ts = 1 / Fs;
    
    fft_single(data,Fs);

    [s, f, t] = spectrogram(data, g, L, Ndft, Fs, 'onesided', 'yaxis');
    s = s ./ M;
    s(2:end-1, :) = 2 * s(2:end-1, :);

    if sc_flag == 1
        load lemi_tf_dat.mat lemi_tf
        fd = (1 / (2 * pi)) * (2 / Ts) * atan((2 * pi * f * Ts) / 2);
        gain = interp1(lemi_tf(:, 1), lemi_tf(:, 2), fd, 'linear', 'extrap');
        phase = interp1(lemi_tf(:, 1), lemi_tf(:, 3), fd, 'linear', 'extrap');
        TF = (gain .* exp(1i * deg2rad(phase)))';
        sc = s ./ (ones(length(t), 1) * TF)';
    else
        sc = s .* 35;
    end

    output.s = sc;
    output.f = f;
    output.t = t;
end


function plotspectrogram(input,cutoff)
    [x,y] = meshgrid(input.t,input.f(input.f<cutoff));
    z = abs(input.s(input.f<cutoff,:));
    z(z>50)=50;
    h = surf(x,y,z);
    get(h);
    set(h,'linestyle','none');
    view(0,90)
    ylim([0,cutoff])
    ylabel('Frequency [Hz]')
    xlim([0,input.t(end)])
    xlabel('Time [s]')
    c = colorbar;
    c.Label.String = ('Magnetic Flux Density [nT]');
    c.Label.Interpreter = 'Latex';
end

function filtered_data = highpass_filter(data, Fs, cutoff_freq)
    [b, a] = butter(4, cutoff_freq / (Fs / 2), 'high');  % 4th order Butterworth
    filtered_data = filtfilt(b, a, data);               % zero-phase filtering
end

function filtered_data = comb_notch_filter(data, Fs, base_freq, num_harmonics, bandwidth)
    filtered_data = data;
    for k = 1:num_harmonics
        notch_freq = k * base_freq;
        Wn = [notch_freq - bandwidth/2, notch_freq + bandwidth/2] / (Fs/2);
        [b, a] = butter(1, Wn, 'stop');
        filtered_data = filtfilt(b, a, filtered_data);
    end
end

function filtered_data = notch_50Hz(input, Fs)
    % Design a comb notch filter for 50 Hz and harmonics - this uses a
    % legacy function
    L = round(Fs / 50);           % spacing between notches
    BW = 1;                       % Bandwidth of each notch (Hz)
    GBW = -5;                     % Gain at the bandwidth edges (dB)
    Nsh = floor(Fs/(2*50)) - 1;   % Number of harmonics to suppress

    d = fdesign.comb('notch','L,BW,GBW,Nsh', L, BW, GBW, Nsh, Fs);
    h = design(d,'SystemObject',true);

    % Apply zero-phase filtering
    filtered_data = filtfilt(h.Numerator, h.Denominator, real(input));
end

function filtered_data = notch_50Hz_iirnotch(data, Fs)
    filtered_data = data;
    base_freq = 50;
    num_harmonics = floor(Fs/(2*base_freq)); % Nyquist limited
    Q = 35; % Quality factor - adjust to control bandwidth
    
    for k = 1:num_harmonics
        f0 = k * base_freq;
        wo = f0 / (Fs/2);  % Normalized frequency
        bw = wo / Q;       % Bandwidth
        [b,a] = iirnotch(wo, bw);
        filtered_data = filtfilt(b, a, filtered_data);
    end
end
function output = fft_single(data,Fs)

    L = length(data);
    f = Fs*(0:(L/2))/L;
    Y = fft(data);
    %matlab fft doesn't normalise by size
    %make double sided spectrum
    P2 = abs(Y/L);
    %make single sided spectrum
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);

    output.data = data;
    output.Y = Y;
    output.P1 = P1;
    output.f = f;
end

function Hd = notch_50_100(Fs, BW)
    % Design notch filter at 50 Hz
    W0_50 = 50 / (Fs/2);
    Q_50 = 50 / BW;
    [b50, a50] = iirnotch(W0_50, W0_50/Q_50);

    % Design notch filter at 100 Hz
    W0_100 = 100 / (Fs/2);
    Q_100 = 100 / BW;
    [b100, a100] = iirnotch(W0_100, W0_100/Q_100);

    % Combine both filters in cascade
    Hd = dsp.FilterCascade( ...
        dsp.IIRFilter('Numerator', b50, 'Denominator', a50), ...
        dsp.IIRFilter('Numerator', b100, 'Denominator', a100));
end