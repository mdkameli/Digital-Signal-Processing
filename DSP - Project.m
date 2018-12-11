clc; clear all; close all;
%STEP 1:
%read the file containing the input signal
DSPproject_054 = 'signal_054.wav';
[x,fs] = audioread(DSPproject_054);

%STEP 2:
%compute the spectrum of the input signal
n = length(x);          % transform length
y = fft(x,n);           % DFT of signal

%STEP 3,4:
%find the frequencies f1 and f2 of the sinusoidal carriers
%find an estimate of the amplitudes A1 and A2 of the carriers;
f = (0:n-1)*(fs/n);
t = (0:1/fs:(length(x)-1)/fs);

Amp = sort(abs(y),'descend');
A1 = 2/n*Amp(1);
A2 = 2/n*Amp(3);
n1 = find(abs(y) == Amp(1));
n2 = find(abs(y) == Amp(3));
f1 = n1(1)*fs/n;
f2 = n2(1)*fs/n;

%STEP 5:
%extracts the two carriers by filtering the input signal
% w1 = ((f1-5)/fs*2);
% w2 = ((f1+5)/fs*2);
% 
% w3 = ((f2-5)/fs*2);
% w4 = ((f2+5)/fs*2);
teta1 = (f2/fs)*2*pi;
delta = (5*pi)/fs;
r = 1-delta;
a1 = -2*r*cos(teta1);
a2 = (r^2);
b0 = delta*2*sin(teta1);
b2 = delta*2*sin(-teta1);
b = [b0 0 b2];
a = [1 a1 a2];

% [b,a]=ellip(2,1,60,[w1,w2]);      % Bandpass digital filter design         
c1 = filter(b,a,x);               % Visualize filter
%%%%
H1 = dsp.IIRFilter('Numerator',b,'Denominator',a);
freqz(H1);
%plot(t,c1);

teta2 = (f1/fs)*2*pi;
delta = (5*pi)/fs;
r = 1-delta;
e1 = -2*r*cos(teta2);
e2 = (r^2);
d0 = delta*2*sin(teta2);
d2 = delta*2*sin(-teta2);
d = [d0 0 d2];
c = [1 e1 e2];
% [d,c]=ellip(2,1,60,[w3,w4]);      % Bandpass digital filter design         
c2 = filter(d,c,x);               % Visualize filter
%%%%
H2 = dsp.IIRFilter('Numerator',d,'Denominator',c);
freqz(H2);
%plot(t,c2);

T1 = 50*c1.*x;
T2 = 50*c2.*x;
%plot(t,T1);

%STEP 7:
%design the filters so that the demodulated signals do not present audible distortions

%LOW-PASS FILTER
Fpass = 4000;            % Passband Frequency
Fstop = 4200;            % Stopband Frequency
Dpass = 0.01;            % Passband Ripple
Dstop = 0.001;           % Stopband Attenuation
% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(fs/2), [1 0], [Dpass, Dstop]);
% Calculate the coefficients using the FIRPM function.
q  = firpm(N, Fo, Ao, W);
Q1 = filter(q,1,T1);
Q2 = filter(q,1,T2);
%%%%
Hd = dfilt.dffir(q);
freqz (Hd);
%plot (t,Q1)

%HIGH-PATH NOTCH
tetan = 0;
rn = 1-((10/fs)*2*pi);
bn1 = -2*cos(tetan);
an1 = -2*rn*cos(tetan);
an2 = (rn^2);
bn = [1 bn1 1];
an = [1 an1 an2];

%%%%
Hn = dsp.IIRFilter('Numerator',bn,'Denominator',an);
freqz (Hn);

y1 = filter(bn,an,Q1);
y2 = filter(bn,an,Q2);

output = [y1 y2];

sound(output, fs);
filename = 'Kameli_Mohammad.wav';
audiowrite(filename,output,fs);

%sound(y1,fs);
%sound(y2,fs);

