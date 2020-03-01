%%
%Drafting File

%Part I.
%Code Given to load sound, visualize, and play
load handel
y = downsample(y,2);
v = y';
vt = fft(v);

%p8 = audioplayer(v,Fs); playblocking(p8);
L=9; n=length(y); 
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k = (2*pi/2*L)*[0:(n/2) -n/2:-1];
ks = fftshift(k);

figure()
subplot(4,1,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

%Gabor Filter
a = 10;
%g = exp(-a*(t-tau).^2);

t_slide = 0:0.1:L;

for j = 1:length(t_slide)
    g = exp(-a*(t-t_slide(j)).^2);
    sg = g.*v; %filtered signal
    sgt = fft(sg); %FFT of filtered signal
    ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal
end


subplot(4,1,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

subplot(4,1,3)
plot(ks, abs(fftshift(sgt)));
xlabel('Frequency'); ylabel('Amplitude');
title('Filtered Signal a = 10');

subplot(4,1,4)
pcolor(t_slide,ks/(2*pi),ssgt.'), 
xlabel('Time'); ylabel('Frequency');
title('Spectrogram of Filtered Signal, a = 10');

shading interp 
%set(gca,'Ylim',[0 50],'Fontsize',16) 
colormap(jet)

%%
%Part I. Multiple Widths
%Code Given to load sound, visualize, and play
load handel
y = downsample(y,2);
v = y';
vt = fft(v);

%p8 = audioplayer(v,Fs); playblocking(p8);
L=9; n=length(y); 
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k = (2*pi/2*L)*[0:(n/2) -n/2:-1];
ks = fftshift(k);


t_slide = 0:0.1:L;

figure()
subplot(2,1,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

subplot(2,1,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

for a = 0:20:100
    figure()

    for j = 1:length(t_slide)
        g = exp(-a*(t-t_slide(j)).^2);
        sg = g.*v; %filtered signal
        sgt = fft(sg); %FFT of filtered signal
        ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal
    end


    subplot(2,1,1)
    plot(ks, abs(fftshift(sgt)));
    xlabel('Frequency'); ylabel('Amplitude');
    title(['Filtered Signal a = ' num2str(a)]);

    subplot(2,1,2)
    pcolor(t_slide,ks/(2*pi),ssgt.'), 
    xlabel('Time'); ylabel('Frequency');
    title(['Spectrogram of Filtered Signal, a = ' num2str(a)]);

    shading interp 
    colormap(jet)
end

%%
%Part I. Changing Sliding
%Code Given to load sound, visualize, and play
load handel
y = downsample(y,2);
v = y';
vt = fft(v);

%p8 = audioplayer(v,Fs); playblocking(p8);
L=9; n=length(y); 
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k = (2*pi/2*L)*[0:(n/2) -n/2:-1];
ks = fftshift(k);


a = 50;

figure()
subplot(2,1,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

subplot(2,1,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

for t_slide_w = 0.1:0.1:1.5
    clear ssgt;
    figure()
    t_slide = 0:t_slide_w:L;
    
    for j = 1:length(t_slide)
        g = exp(-a*(t-t_slide(j)).^2);
        sg = g.*v; %filtered signal
        sgt = fft(sg); %FFT of filtered signal
        ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal
    end


    subplot(2,1,1)
    plot(ks, abs(fftshift(sgt)));
    xlabel('Frequency'); ylabel('Amplitude');
    title(['Filtered Signal Time Width = ' num2str(t_slide_w)]);

    subplot(2,1,2)
    pcolor(t_slide,ks/(2*pi),ssgt.'), 
    xlabel('Time'); ylabel('Frequency');
    title(['Spectrogram of Filtered Signal, Time Width = ' num2str(t_slide_w)]);

    shading interp 
    colormap(jet)
end

%%
%Mexican Hat Filter
load handel
y = downsample(y,2);
v = y';
vt = fft(v);

%p8 = audioplayer(v,Fs); playblocking(p8);
L=9; n=length(y); 
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k = (2*pi/2*L)*[0:(n/2) -n/2:-1];
ks = fftshift(k);

% tau = 4;
% a = 50;
% m_hat = (1-a*(t-tau).^2).*exp((-1*a*(t-tau).^2)/2);
% figure()
% plot(m_hat)

t_slide = 0:0.1:L;
a = 1000;

for j = 1:length(t_slide)
    m_hat = (1-a*(t-t_slide(j)).^2).*exp((-1*a*(t-t_slide(j)).^2)/2);
    sg = m_hat.*v; %filtered signal
    sgt = fft(sg); %FFT of filtered signal
    ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal
end

figure()
subplot(4,1,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

subplot(4,1,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

subplot(4,1,3)
plot(ks, abs(fftshift(sgt)));
xlabel('Frequency'); ylabel('Amplitude');
title('Signal Filtered with Mexican Hat Wavelet');

subplot(4,1,4)
pcolor(t_slide,ks/(2*pi),ssgt.'), 
xlabel('Time'); ylabel('Frequency');
title('Spectrogram of Mexican Hat Filtered Signal');
ylim ([0 60000]);

shading interp  
colormap(jet)

%%
%Shannon Window


load handel
y = downsample(y,2);
v = y';
vt = fft(v);

%p8 = audioplayer(v,Fs); playblocking(p8);
L=9; n=length(y); 
t2=linspace(0,L,n+1); 
t=t2(1:n); 
k = (2*pi/2*L)*[0:(n/2) -n/2:-1];
ks = fftshift(k);


%m_hat = (1-(t-t_slide)^2)*exp((-1*(t-t_slide)^2)/2);

t_slide = 0:0.1:L;

tau = 4;
width = 0.1;

%shannon = heaviside(t-(tau - width/2)) - heaviside(t-(tau + width/2));

for j = 1:length(t_slide)
    sg = (heaviside(t-(t_slide(j) - width/2)) - heaviside(t-(t_slide(j) + width/2))).*v; %filtered signal
    sgt = fft(sg); %FFT of filtered signal
    ssgt(j,:) = abs(fftshift(sgt)); %Shifted filter Signal
end

figure()
subplot(4,1,1)
plot(t,v);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Signal of Interest , v(n)');

subplot(4,1,2)
plot(ks, abs(fftshift(vt)));
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Unfiltered Signal');

subplot(4,1,3)
plot(ks, abs(fftshift(sgt)));
xlabel('Frequency'); ylabel('Amplitude');
title('Signal Filtered with Shannon Window');

subplot(4,1,4)
pcolor(t_slide,ks/(2*pi),ssgt.'), 
xlabel('Time'); ylabel('Frequency');
title('Spectrogram of Shannon Window Filtered Signal');
%ylim([0 30]);

shading interp  
colormap(jet)

%%
%Plotting Functions

tau = 4;
a = 50;
width = 3;

G = exp(-a*(t-tau).^2);
m_hat = (1-a*(t-tau).^2).*exp((-1*a*(t-tau).^2)/2);
shannon = heaviside(t-(tau - width/2)) - heaviside(t-(tau + width/2));

figure()
subplot(3,1,1)
plot(t,G)
axis([-inf inf -0.5 1.1])
title('Gaussian Filter')

subplot(3,1,2)
plot(t,m_hat)
axis([-inf inf -0.5 1.1])
title('Mexican Hat Wavelet')

subplot(3,1,3)
plot(t,shannon)
axis([-inf inf -0.5 1.2])
title('Shannon Window')


