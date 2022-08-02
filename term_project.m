%Design and Implement instantaneous amplitude envelope estimator 
%by generating analytical signal using the Hilbert transformer. 
%Use DFT for design of Hilbert transformer. Use heart sound signal 
%with sampling rate of 2000 Hz. and Assume block duration of 10 second

fs = 2e3;                                   % Sampling Frequency in Hz

[y, Fs] = audioread("beats.wav");           % Obtaining the signal from the audio file
y = y';

t = (0:length(y)-1)/Fs;                     % Time Domain of the obtained signal

y_fft = R2_DIT(y, -1);                      % Obtaining the analytical signal
y_fft_ideal = fft(y);                       % Ideal FFT for verification

complex_y_fft = 1j*y_fft;
n = length(y_fft);

positive_frequencies = 2:floor(n/2)+mod(n, 2);
negative_frequencies = ceil(n/2)+1+~mod(n, 2):n;

y_fft(positive_frequencies) = y_fft(positive_frequencies) + -1j*complex_y_fft(positive_frequencies);
y_fft(negative_frequencies) = y_fft(negative_frequencies) +  1j*complex_y_fft(negative_frequencies);

z = R2_DIT(y_fft, 1);

inst_amp = abs(z);                          % Obtaining instantaneous amplitude envelope

z_ideal = hilbert(y);                       % Obtaining the analytical signal using in-built function

inst_amp_ideal = abs(z_ideal);              % Obtaining instantaneous amplitude envelope of ideal analytical signal

L = 100;                                    % Designing an MA filter to smoothen the envelope
b = (1/L)*ones(1, L);
a = 1;

inst_amp_smooth = zeros(n, 1);

for i = 1:n
    temp = 0;
    for j = 1:L
        if (i-j > 0)
            temp = temp + inst_amp(i-j);
        end
    end
    inst_amp_smooth(i) = temp/L;
    temp = 0;
end

delay = (L-1)/2;                            % Delay introduced due to the MA filter

a0 = 2*fs;                                  % Taking a section of the entire signal for better legibility
b0 = 4*fs;

figure
plot(t(a0:b0), y(a0:b0));
hold on;
plot(t(a0:b0)-(delay/fs), inst_amp_smooth(a0:b0), "r");
xlabel("Time (s)");
xlim([a0/fs b0/fs]);
title("Smoothened Instantaneous Amplitude Envelope using Moving Average Filter");

figure
plot(t(a0:b0), y(a0:b0));
hold on;
plot(t(a0:b0), inst_amp(a0:b0), "r");
xlabel("Time (s)");
title("Instantaneous Amplitude Envelope using Hilbert's Transform");

figure
plot(t(a0:b0), y(a0:b0));
hold on;
plot(t(a0:b0), inst_amp(a0:b0), "g");
hold on;
plot(t(a0:b0), inst_amp_ideal(a0:b0), "r");
title("Instantaneous Amplitude Envelope from both methods");
xlabel("Time (s)");
legend("Sampled Signal", "Written Code", "In-built Function", "location", "best");

figure
plot(t(a0:b0), y(a0:b0));
hold on;
plot(t(a0:b0), inst_amp_ideal(a0:b0), "r");
xlabel("Time (s)");
title("Instantaneous Amplitude Envelope using in-built function");

figure
plot(t(a0:b0), y(a0:b0));
xlabel("Time (s)");
title("The Section of Interest");

figure
plot(t, y);
xlabel("Time (s)");
title("Heart Sound Signal sampled at 2kHz");

function [y] = R2_DIT(x, type)              %function computes FFT when type = -1 and IFFT when type = 1
p=2^nextpow2(length(x));
x=[x zeros(1,p-length(x))];
N=length(x);
S=log2(N);
Half=1;
x=bitrevorder(x);
for stage=1:S
    for index=0:(2^stage):(N-1)
        for n=0:(Half-1)
            pos=n+index+1;
            pow=(2^(S-stage))*n;
            w=exp((type*1j)*(2*pi)*pow/N);
            a=x(pos)+x(pos+Half).*w;
            b=x(pos)-x(pos+Half).*w;
            x(pos)=a;
            x(pos+Half)=b;
        end
    end
Half=2*Half;
end
if type == 1
    y = x/N;
else
    y = x;
end
end

function [y] = DFT(x, type)      %#ok<*DEFNU> %function computes DFT when type = -1 and IDFT when type = 1
    N = length(x);
    y = zeros(N, 1);
    for k = 1:N-1
        for n = 1:N
            y(k+1) = y(k)+x(n)*exp(type*1j*2*pi*n*k/N);
        end
    end
    if type == 1
       y = y/N; 
    end
end