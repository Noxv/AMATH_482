%Max Walter
%AMATH 482
%Homework 2: Gabor Transform and Signal Processing
%2/7/2020
%University of Washington
clear all; close all; clc

%Part I. Looking at different filters manipulating width and time window.

%Loading in built in music file of Handel's Hallelujah Chorus

load handel %Outputs [y Fs]
y = downsample(y,2); %Downsampling for improved performance
v = y'; %Comes out as a column vector
vt = fft(v); %Initial Fourier Transform of the Data

%p8 = audioplayer(v,Fs); playblocking(p8); %Code uses to play the file

%Defining Parameters
L=9; n=length(y); 
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k = (1/L)*[0:(n/2) -n/2:-1];
ks = fftshift(k);

%Plotting Original Time Signal
figure()
subplot(2,2,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

%Plotting Unfiltered Signal
subplot(2,2,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

%Gaussian Filter
%g = exp(-a*(t-tau).^2);

a = 10; %Width
t_slide = 0:0.1:L; %Sliding Window

for j = 1:length(t_slide)
    g = exp(-a*(t-t_slide(j)).^2);
    sg = g.*v; %Filtered Signal
    sgt = fft(sg); %FFT of Filtered Signal
    ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal for making spectrogram.
end

%Plotting Filtered Signal
subplot(2,2,3)
plot(ks, abs(fftshift(sgt)));
xlabel('Frequency'); ylabel('Amplitude');
title('Filtered Signal, a = 10');

%Plotting Spectrogram
subplot(2,2,4)
pcolor(t_slide,ks,ssgt.'), 
xlabel('Time'); ylabel('Frequency');
title('Spectrogram of Filtered Signal, a = 10');
set(gca,'Ylim',[0 Fs/4]) 

shading interp 
colormap(jet)

%%
%%%%%
%Plotting Different Width


%Plotting Original Signal
figure()
subplot(2,1,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

%Plotting Unfiltered
subplot(2,1,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

%Running through different Widths
IDX = 1;
figure()
for a = 0:25:100

    for j = 1:length(t_slide)
        g = exp(-a*(t-t_slide(j)).^2);
        sg = g.*v;
        sgt = fft(sg);
        ssgt(j,:) = abs(fftshift(sgt)); 
    end

    subplot(2,5,IDX)
    plot(ks, abs(fftshift(sgt)));
    xlabel('Frequency'); ylabel('Amplitude');
    title(['Filtered Signal a = ' num2str(a)]);
    
    IDX = IDX + 5;
    subplot(2,5,IDX)
    pcolor(t_slide,ks,ssgt.'), 
    xlabel('Time'); ylabel('Frequency');
    title(['Spectrogram of Filtered Signal, a = ' num2str(a)]);
    
    IDX = IDX - 4;
    shading interp
    set(gca,'Ylim',[0 Fs/4]) 
    colormap(jet)
end

%%
%%%%%%%%
%Testing Different Sized Sampling Windows
%Looking at sizes from 0.1 second to 1.5 seconds

a = 50; %Redefining filterd width from previous section
IDX = 1;
figure()
%Filtering Times
for t_slide_w = 0.1:0.2:1.1
    clear ssgt;

    t_slide = 0:t_slide_w:L;
    
    %Running New slide over signal through Gaussian Filter
    for j = 1:length(t_slide)
        g = exp(-a*(t-t_slide(j)).^2);
        sg = g.*v; %filtered signal
        sgt = fft(sg); %FFT of filtered signal
        ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal
    end

    %Plotting Filtered Signal
    subplot(2,6,IDX)
    plot(ks, abs(fftshift(sgt)));
    xlabel('Frequency'); ylabel('Amplitude');
    title(['Filtered Signal Time Width = ' num2str(t_slide_w)]);
    IDX = IDX + 6;

    %Plotting Unfiltered Signal
    subplot(2,6,IDX)
    pcolor(t_slide,ks,ssgt.'), 
    xlabel('Time'); ylabel('Frequency');
    title(['Spectrogram of Filtered Signal, Time Width = ' num2str(t_slide_w)]);
    IDX = IDX - 5;
    
    shading interp
    set(gca,'Ylim',[0 Fs/4])  
    colormap(jet)

end

%%
%%%%%%
%Mexican Hat Filter

%Using Values from previous tests
t_slide = 0:0.1:L;
a = 50;


for j = 1:length(t_slide)
    m_hat = (1-a*(t-t_slide(j)).^2).*exp((-1*a*(t-t_slide(j)).^2)/2); %Mexican Hat Function
    sg = m_hat.*v;
    sgt = fft(sg);
    ssgt(j,:) = abs(fftshift(sgt));
end


figure()
%Filtered Signal with Mexican Hat Wavelet
subplot(1,2,1)
plot(ks, abs(fftshift(sgt)));
xlabel('Frequency'); ylabel('Amplitude');
title('Signal Filtered with Mexican Hat Wavelet');

subplot(1,2,2)
pcolor(t_slide,ks,ssgt.'), 
xlabel('Time'); ylabel('Frequency');
title('Spectrogram of Mexican Hat Filtered Signal');

shading interp
set(gca,'Ylim',[0 Fs/4])  
colormap(jet)


%%
%%%%%%%%%
%Shannon Window
%Heaviside Function to Step and Act as Shannon window
%First part ramps up holds at 1, then ramps down

t_slide = 0:0.1:L;
width = 0.1; %Selecting hold width of step function

for j = 1:length(t_slide)
    sg = (heaviside(t-(t_slide(j) - width/2)) - heaviside(t-(t_slide(j) + width/2))).*v; 
    sgt = fft(sg); %FFT of filtered signal
    ssgt(j,:) = abs(fftshift(sgt));
end

figure()
%Signal after Shannon Window
subplot(1,2,1)
plot(ks, abs(fftshift(sgt)));
xlabel('Frequency'); ylabel('Amplitude');
title('Signal Filtered with Shannon Window');

%Spectrogram of Signal
subplot(1,2,2)
pcolor(t_slide,ks,ssgt.'), 
xlabel('Time'); ylabel('Frequency');
title('Spectrogram of Shannon Window Filtered Signal');
set(gca,'Ylim',[0 Fs/4])  

shading interp  
colormap(jet)

%%
%%%%
%Plotting Commands for functions used to filter

tau = 4;
a = 50;
width = 3;

%Gaussian
G = exp(-a*(t-tau).^2);

%Mexican Hat
m_hat = (1-a*(t-tau).^2).*exp((-1*a*(t-tau).^2)/2);

%Shannon
shannon = heaviside(t-(tau - width/2)) - heaviside(t-(tau + width/2));

figure()
subplot(3,2,1)
plot(t,G,'Linewidth',3)
axis([-inf inf -0.5 1.1])
title('Gaussian Filter')

subplot(3,2,2)
plot(t,G.*v)
hold on
plot(t,G,'Linewidth',3)
axis([-inf inf -0.5 1.1])
title('Signal Filtered by Gaussian Filter')

subplot(3,2,3)
plot(t,m_hat,'Linewidth',3)
axis([-inf inf -0.5 1.1])
title('Mexican Hat Wavelet')

subplot(3,2,4)
plot(t,m_hat.*v)
hold on
plot(t,m_hat,'Linewidth',3)
axis([-inf inf -0.5 1.1])
title('Signal Filtered by Mexican Hat Wavelet')

subplot(3,2,5)
plot(t,shannon,'Linewidth',3)
axis([-inf inf -0.5 1.2])
title('Shannon Window')

subplot(3,2,6)
plot(t,shannon.*v)
hold on
plot(t,shannon,'Linewidth',3)
axis([-inf inf -0.5 1.2])
title('Shannon Window')


