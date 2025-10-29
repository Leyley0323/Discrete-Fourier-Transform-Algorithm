% DSP Discrete Fourier Transform Project using MATLAB

clear; close all; clc;

% function to implement DFT
function [X, f] = myDFT(x, Fs)
    N = length(x); % gets # of samples in x
    X = zeros(1, N); % for zero padding & starting X at 0

    % DFT Forumla
    for k = 0:N-1
        for n = 0:N-1
            X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
        end
    end

    % returns positive frequencies from 0 to Fs/2
    f = (0:(N/2)-1) * (Fs / N);
    X = X(1:N/2);
end

%case (a): 1 cycle of x(t)=47cos(2pi200t), 20 zeros, Fs=2kHz
Fs1 = 2000; 
Ts1 = 1/Fs1; %sampling period
fc1 = 200;
tc1 = 1/fc1; 
Nsamp = round(Fs1/fc1); % Samples per cycle
ta = (0:Nsamp-1) * Ts1;
signal1 = 47*cos(2*pi*fc1*ta);
xa = [signal1, zeros(1, 20)]; % 20 zeros padded

[Xa, fa] = myDFT(xa, Fs1);
AmpA = (2 * abs(Xa)) / Nsamp;


figure(1);
plot(fa, AmpA, 'c-', 'LineWidth', 2);
title('Case A: 1 cycle of x(t) = 47cos(2pi200t) with 20 zeros');
xlabel('Frequency');
ylabel('Amplitude');
grid on;
xlim([0 Fs1/2]); % Nyquist sampling rate

%case (b): 1 cycle of x(t)=47cos(2pi200t), 200 zeros, Fs=2kHz 
Fs2 = 2000; 
Ts2 = 1/Fs2; %sampling period
fc2 = 200;
tc2 = 1/fc2; 
Nsamp2 = round(Fs2/fc2); % Samples per cycle
tb = (0:Nsamp2-1) * Ts2;

signal2 = 47*cos(2*pi*fc2*tb);
xb = [signal2, zeros(1, 200)]; % 200 zeros padded

[Xb, fb] = myDFT(xb, Fs2);
len2 = length(signal2); 
AmpB = (2 * abs(Xb)) / Nsamp2;

figure(2);
plot(fb, AmpB, 'm-', 'LineWidth', 2);
title('Case B: 1 cycle of x(t) = 47cos(2pi200t) with 200 zeros');
xlabel('Frequency');
ylabel('Amplitude');
grid on;
xlim([0 Fs2/2]);

%case (c): 100 cycles of x(t)=47cos(2pi200t), 200 zeros, Fs=2kHz
Fs3 = 2000; 
Ts3 = 1/Fs3; %sampling period
fc3 = 200;
tc3 = 100/fc3; %100 cycles
tc = 0:Ts3:tc3-Ts3;
signal3 = 47*cos(2*pi*fc3*tc);
xc = [signal3, zeros(1, 200)]; % 200 zeros padded

[Xc, fc] = myDFT(xc, Fs3);
len3 = length(signal3); 
AmpC = (2 * abs(Xc)) / len3;

figure(3);
plot(fc, AmpC, 'c-', 'LineWidth', 2);
title('Case C: 100 cycles of x(t) = 47cos(2pi200t) with 200 zeros');
xlabel('Frequency');
ylabel('Amplitude');
grid on;
xlim([0 Fs3/2]);

%case (d) x(t)=50cos(2pi200t)sin(2pi200t), 4000 zeros, 
% Fs= 2kHz, 0<= t <= 1s
Fs4 = 2000; 
Ts4 = 1/Fs4; %sampling period
fc4 = 200; 
td = 0:Ts4:1-Ts4; % 1 - (1/2k) = 0.9995 seconds 
signal4 = 50*cos(2*pi*fc4*td) .* sin(2*pi*fc4*td);
xd = [signal4, zeros(1, 4000)]; % 4000 zeros padded

[Xd, fd] = myDFT(xd, Fs4);
len4 = length(signal4); 
AmpD = (2 * abs(Xd)) / len4;

figure(4);
plot(fd, AmpD, 'm-', 'LineWidth', 2);
title('Case D: x(t) = 50cos(2pi200t)sin(2pi200t) with 4000 zeros');
xlabel('Frequency');
ylabel('Amplitude');
grid on;
xlim([0 Fs4/2]);

%case (e) x(t)=47cos(2pi300t) + 50sin(2pi400t) + 74cos(2pi1500t) +
% 5cos(2pi1600t), 8000 zeros, Fs= 4kHz, 0 <= t <= 0.5s
Fs5 = 4000; 
Ts5 = 1/Fs5; %sampling period
te = 0:Ts5:0.5-Ts5; % (1/2) - (1/4k) = 0.49975 seconds 
signal5 = 47*cos(2*pi*300*te) + 50*sin(2*pi*400*te) + 74*cos(2*pi*1500*te) + 5*cos(2*pi*1600*te);
xe = [signal5, zeros(1, 8000)]; % 8000 zeros padded

[Xe, fe] = myDFT(xe, Fs5);
len5 = length(signal5); 
AmpE = (2 * abs(Xe)) / len5;

figure(5);
plot(fe, AmpE, 'c-', 'LineWidth', 2);
title('Case E: x(t)=47cos(2pi300t) + 50sin(2pi400t) + 74cos(2pi1500t) + 5cos(2pi1600t), Fs = 4kHz');
xlabel('Frequency');
ylabel('Amplitude');
grid on;
xlim([0 Fs5/2]);

%case (f) x(t)=47cos(2pi300t) + 50sin(2pi400t) + 74cos(2pi1500t) +
% 5cos(2pi1600t), 8000 zeros, Fs= 1.2kHz, 0 <= t <= 0.5s
Fs6 = 1200; 
Ts6 = 1/Fs6; %sampling period
tf = 0:Ts6:0.5-Ts6; % (1/2) - (1/1.2k) = 0.49917 seconds 
signal6 = 47*cos(2*pi*300*tf) + 50*sin(2*pi*400*tf) + 74*cos(2*pi*1500*tf) + 5*cos(2*pi*1600*tf);
xf = [signal6, zeros(1, 8000)]; % 8000 zeros padded

[Xf, ff] = myDFT(xf, Fs6);
len6 = length(signal6); 
AmpF = (2 * abs(Xf)) / len6;

figure(6);
plot(ff, AmpF, 'm-', 'LineWidth', 2);
title('Case F: x(t)=47cos(2pi300t) + 50sin(2pi400t) + 74cos(2pi1500t) + 5cos(2pi1600t), Fs = 1.2kHz');
xlabel('Frequency');
ylabel('Amplitude');
grid on;
xlim([0 Fs6/2]);

