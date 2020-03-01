%Max Walter
%AMATh 482
%Homework 1
%1/15/2020

clear all; close all; clc;

load('Testdata.mat')


L=15; % spatial domain
n=64; % Fourier modes

%Space
x2=linspace(-L,L,n+1); 
x=x2(1:n); 
y=x; 
z=x; 
[X,Y,Z]=meshgrid(x,y,z); %64x64x64 Grid for x,y,z spatial location

%Frequency
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
ks=fftshift(k);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); %Grid of Frequency from -2pi to 2pi

%Using Avg. to cancel/reduce whitenoise
fft_sum = zeros(64,64,64);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    fft_run = fftn(Un);
    fft_sum = fft_sum + fft_run;
end
abs_fft_sum = abs(fftshift(fft_sum/20));

%Finding maximum value and index in avg. data
[value,index] = max(abs_fft_sum(:));

%using index to find subscripts
[sX,sY,sZ] = ind2sub([n,n,n],index);

%Converting to frequency
Fx = Kx(sX,sY,sZ);
Fy = Ky(sX,sY,sZ);
Fz = Kz(sX,sY,sZ);

%Creating Filter using Gaussian Filter from book
filter = exp(-0.2 * ((Kx - Fx).^2 + (Ky - Fy).^2 + (Kz - Fz).^2));

%Initalizing Positions
x_pos = zeros(1,20); 
y_pos = zeros(1,20); 
z_pos = zeros(1,20);

%Finding Marble Position
for j=1:20
    
    %Filtering Signal
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    fft_run = fftn(Un);
    filtered_signal = filter .* fftshift(fft_run);
    
    %Returning to Spacial
    Clean_Position = real(ifftn(filtered_signal));
    
    %Finding Marble using largest signal like in previous method
    [value,index] = max(Clean_Position(:));
    [sX,sY,sZ] = ind2sub([n,n,n],index);
    
    %Storing Location
    x_pos(j) = X(sX,sY,sZ);
    y_pos(j) = Y(sX,sY,sZ);
    z_pos(j) = Z(sX,sY,sZ);
end

%Ultra sonic strike
endpoint = [x_pos(end), y_pos(end), z_pos(end)];

figure(1)
plot3(x_pos,y_pos,z_pos, 'Linewidth', 5)
hold on
plot3(x_pos(end), y_pos(end), z_pos(end),'ro', 'Linewidth', 10)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
legend('Path','Ultrasonic Strike')
title('Path of Marble in Fluffy')
