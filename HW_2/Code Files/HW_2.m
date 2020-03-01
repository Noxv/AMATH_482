%AMATH 482 HW 2
%Code by Max Walter
%1/30/2020 V1 Draft


%%
%Part II.
clear all; close all; clc

%Piano
%figure()
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds 
L = round(length(y)/Fs);
y = downsample(y,5);

% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude'); 
% title('Mary had a little lamb (piano)');
%p8 = audioplayer(y,Fs); playblocking(p8);

piano_signal = y';
n = length(y);
t2p = linspace(0,tr_piano,n+1);
tp=t2p(1:n); 

kp = (2*pi/L)*[0:(n/2)-1 -n/2:-1];
ksp = fftshift(kp);

t_slide_p = 0:0.1:L;
a = 100;


pnotes = zeros(length(t_slide_p),1); %Initilizing pnotes vector
spp = zeros(length(t_slide_p),length(y)); %Initilizing Shifted Signal Matrix

for j = 1:length(t_slide_p)
    g = exp(-a*(tp - t_slide_p(j)).^2);
    fps = g.*piano_signal; %filtered piano signal
    lfps = lowpass(fps,370,Fs);
    fpst = fft(lfps); %filtered Piano Signal Transform
    
    [M,I] = max(fpst); %Finding Index where signal is maximum for note
    pnotes(j,:) = abs(kp(I))/(2*pi); %Using Finding maximum index in frequency.
    spp(j,:) = abs(fftshift(fpst)); %Shifting signal for spectrogram
end

% figure()
% plot(ksp, fftshift(fpst))

figure();
subplot(2,2,1)
pcolor(t_slide_p , ksp/(2*pi) ,spp.') 
shading interp 
set(gca,'Ylim',[0 1000],'Fontsize',16) 
colormap(hot)
xlabel("Time (s)");
ylabel("Frequency (Hz)"); 
title("Piano");


%Recorder

[y,Fs] = audioread('music2.wav');
tr_recorder=length(y)/Fs; % record time in seconds 
Lr = round(length(y)/Fs);
y = downsample(y,5);

% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude'); 
% title('Mary had a little lamb (piano)');
%p8 = audioplayer(y,Fs); playblocking(p8);

recorder_signal = y';
n = length(y);
t2r = linspace(0,tr_recorder,n+1);
tr=t2r(1:n); 

kr = (2*pi/Lr)*[0:(n/2) -n/2:-1];
ksr = fftshift(kr);

t_slide_r = 0:0.1:Lr;
a = 100;


rnotes = zeros(length(t_slide_r),1);
spr = zeros(length(t_slide_r),length(y));

for j = 1:length(t_slide_r)
    g = exp(-a*(tr - t_slide_r(j)).^2);
    frs = g.*recorder_signal; %filtered recorder signal
    frst = fft(frs); %filtered Recorder Signal Transform
    
    [M,I] = max(frst);
    rnotes(j,:) = abs(kp(I))/(2*pi);
    spr(j,:) = abs(fftshift(frst));
end


subplot(2,2,3)
pcolor(t_slide_r , ksr/(2*pi) ,spr.') 
shading interp 
set(gca,'Ylim',[0 1500],'Fontsize',16) 
colormap(hot)
xlabel("Time (s)");
ylabel("Frequency (Hz)"); 
title("Recorder");


subplot(2,2,2)
plot(t_slide_p,pnotes+10,'o','MarkerFaceColor', 'k'); 
yticks([246.94, 261.63, 277.18, 293.66, 311.13, 329.63, 349.23]); 
yticklabels({'B3','C4','C#4','D4','E4','F4'});
ylim ([246 350])
title("Score for Piano");
xlabel("Time (s)"); ylabel("Note");

subplot(2,2,4)
plot(t_slide_r,rnotes + 20,'o', 'MarkerFaceColor', 'k') 
yticks([698.46, 739.99, 783.99, 830.61, 880, 932.33]); 
yticklabels({'F5','F#5','G5','G#5','A5','A#5'});
ylim ([680 950])
title("Score for Recorder");
xlabel("Time (s)"); ylabel("Note");


