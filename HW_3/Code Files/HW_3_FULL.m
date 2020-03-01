%AMATH 482
%MAX WALTER
%HOMEWORK 3
clear all; close all; clc
tic
1
%Loading all data
load('cam1_1.mat')
load('cam2_1.mat')
load('cam3_1.mat')
numFrames1 = size(vidFrames1_1,4);
numFrames2 = size(vidFrames2_1,4);
numFrames3 = size(vidFrames3_1,4);

%%
%Data 1

%Creating Frame for bucket
frame = zeros(480,640);
frame(200:430,300:400) = 1;
frame_u8 = uint8(frame);


%Initializing In-Loop Values
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));


for j = 1:numFrames1
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X1 = vidFrames1_1(:,:,:,j); %Loading Image Frames
    X1_g = rgb2gray(X1); %Converting to Gray
    fX1_g = X1_g .* frame_u8; %Frame only the bucket moving.

    %Focus on the Light by making B/W
    light1 = fX1_g > 250;
    idx1 = find(light1);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light1),idx1);
    
    %Since there are two points, finding the mean point for the can
    mean_x1(j) = mean(X);
    mean_y1(j) = mean(Y);
end

%Inverting image
corrected_image1 = 480 - mean_y1;

%Plotting
figure()
hold on
plot(1:numFrames1,corrected_image1)

%%
%Data 2

%Creating Frame for bucket
frame = zeros(480,640);
frame(100:385,260:340) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));


for j = 1:numFrames2
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X2 = vidFrames2_1(:,:,:,j); %Loading Image Frames
    X2_g = rgb2gray(X2); %Converting to Gray
    fX2_g = X2_g .* frame_u8; %Frame only the bucket moving.

   
    %Focus on the Light by making B/W
    light2 = fX2_g > 240;
    idx2 = find(light2);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light2),idx2);
    
    %Since there are two points, finding the mean point for the can
    mean_x2(j) = mean(X);
    mean_y2(j) = mean(Y);
end

%Inverting image
corrected_image2 = 480 - mean_y2;

%Plotting
plot(1:numFrames2,corrected_image2)

%%
%Data 3

%Creating Frame for bucket
frame = zeros(480,640);
frame(190:330,270:485) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));


for j = 1:numFrames3
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X3 = vidFrames3_1(:,:,:,j); %Loading Image Frames
    X3_g = rgb2gray(X3); %Converting to Gray
    fX3_g = X3_g .* frame_u8; %Frame only the bucket moving.
    
    %Focus on the Light by making B/W
    light3 = fX3_g > 235;
    idx3 = find(light3);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light3),idx3);
    
    %Since there are two points, finding the mean point for the can
    mean_x3(j) = mean(X);
    mean_y3(j) = mean(Y);
end

%Inverting image, using X component as the video is titled
corrected_image3 = mean_x3;

%Plotting
plot(1:numFrames3,corrected_image3)
title('Uncentered Plots')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')

%%
%All data

%Alligning Data
[M,I1] = min(corrected_image1(1,[1:75]));
min1 = I1;

[M,I2] = min(corrected_image2(1,[1:75]));
min2 = I2;

[M,I3] = min(corrected_image3(1,[1:75]));
min3 = I3;



width = 191;
shift_img1 = 480 - mean_y1(min1:min1+width);
shift_img2 = 480 - mean_y2(min2:min2+width);
shift_img3 = 480 - mean_x3(min3:min3+width);

figure()
hold on
plot(1:width+1,shift_img1)
plot(1:width+1,shift_img2)
plot(1:width+1,shift_img3)
title('Alligned Signals')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')

Full_set = [mean_x1(min1:min1+width); shift_img1; mean_x2(min2:min2+width); shift_img2; shift_img3; mean_y3(min3:min3+width);];
[m,n] = size(Full_set);
mean_FS = mean(Full_set,2);
Full_set = Full_set - repmat(mean_FS,1,n);

% figure()
% plot(1:192,mean_FS)
% title('Mean Signal')
% xlabel('Frames')
% ylabel('Displacement')
% legend('Mean')

[U,S,V] = svd(Full_set/sqrt(n-1));

lambda = diag(S).^2;
sigma = diag(S);
energy = lambda/sum(lambda);

proj = U'*Full_set;

figure()
plot(1:6, [energy(1) sum(energy(1:2)) sum(energy(1:3)) sum(energy(1:4)) sum(energy(1:5)) sum(energy(1:6))],'mo')
title('Primary Component Recreation')
xlabel('Frames')
ylabel('Displacement')
legend('PC1')

figure()
plot(1:192,proj(1,:))
title('Primary Component Recreation')
xlabel('Frames')
ylabel('Displacement')
legend('PC1')

%%
clear all;
2
%Loading all data
load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')
numFrames1 = size(vidFrames1_2,4);
numFrames2 = size(vidFrames2_2,4);
numFrames3 = size(vidFrames3_2,4);

%%
%Data 1

%Creating Frame for bucket
frame = zeros(480,640);
frame(200:430,300:400) = 1;
frame_u8 = uint8(frame);


%Initializing In-Loop Values
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));


for j = 1:numFrames1
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X1 = vidFrames1_2(:,:,:,j); %Loading Image Frames
    X1_g = rgb2gray(X1); %Converting to Gray
    fX1_g = X1_g .* frame_u8; %Frame only the bucket moving.

    %Focus on the Light by making B/W
    light1 = fX1_g > 250;
    idx1 = find(light1);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light1),idx1);
    
    %Since there are two points, finding the mean point for the can
    mean_x1(j) = mean(X);
    mean_y1(j) = mean(Y);
end

%Inverting image
corrected_image1 = 480 - mean_y1;

%Plotting
figure()
hold on
plot(1:numFrames1,corrected_image1)

%%
%Data 2

%Creating Frame for bucket
frame = zeros(480,640);
frame(100:430,100:460) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));


for j = 1:numFrames2
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X2 = vidFrames2_2(:,:,:,j); %Loading Image Frames
    X2_g = rgb2gray(X2); %Converting to Gray
    fX2_g = X2_g .* frame_u8; %Frame only the bucket moving.

   
    %Focus on the Light by making B/W
    light2 = fX2_g > 240;
    idx2 = find(light2);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light2),idx2);
    
    %Since there are two points, finding the mean point for the can
    mean_x2(j) = mean(X);
    mean_y2(j) = mean(Y);
end

%Inverting image
corrected_image2 = 480 - mean_y2;

%Plotting
plot(1:numFrames2,corrected_image2)

%%
%Data 3

%Creating Frame for bucket
frame = zeros(480,640);
frame(190:330,270:485) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));


for j = 1:numFrames3
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X3 = vidFrames3_2(:,:,:,j); %Loading Image Frames
    X3_g = rgb2gray(X3); %Converting to Gray
    fX3_g = X3_g .* frame_u8; %Frame only the bucket moving.
    
    %Focus on the Light by making B/W
    light3 = fX3_g > 235;
    idx3 = find(light3);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light3),idx3);
    
    %Since there are two points, finding the mean point for the can
    mean_x3(j) = mean(X);
    mean_y3(j) = mean(Y);
end

%Inverting image, using X component as the video is titled
corrected_image3 = 480 - mean_x3;

%Plotting
plot(1:numFrames3,corrected_image3)
title('Uncentered Plots')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')

%%
%All data

%Alligning Data
[M,I1] = min(corrected_image1(1,[1:75]));
min1 = I1;

[M,I2] = min(corrected_image2(1,[1:75]));
min2 = I2;

[M,I3] = min(corrected_image3(1,[1:75]));
min3 = I3;



width = 191;
shift_img1 = 480 - mean_y1(min1:min1+width);
shift_img2 = 480 - mean_y2(min2:min2+width);
shift_img3 = 480 - mean_x3(min3:min3+width);

figure()
hold on
plot(1:width+1,shift_img1)
plot(1:width+1,shift_img2)
plot(1:width+1,shift_img3)
title('Alligned Signals')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')


Full_set = [mean_x1(min1:min1+width); shift_img1; mean_x2(min2:min2+width); shift_img2; shift_img3; mean_y3(min3:min3+width)];
[m,n] = size(Full_set);
mean_FS = mean(Full_set,2);
Full_set = Full_set - repmat(mean_FS,1,n);

% figure()
% plot(1:192,mean_FS)
% title('Mean Signal')
% xlabel('Frames')
% ylabel('Displacement')
% legend('Mean')

[U,S,V] = svd(Full_set/sqrt(n-1));

lambda = diag(S).^2;
sigma = diag(S);
energy = lambda/sum(lambda);

proj = U'*Full_set;

figure()
plot(1:6, [energy(1) sum(energy(1:2)) sum(energy(1:3)) sum(energy(1:4)) sum(energy(1:5)) sum(energy(1:6))],'mo')
title('Energy of each Diagonal Variance')
xlabel('Diagonal Variances')
ylabel('Energy Captured')

figure()
plot(1:192,proj(1,:))
title('Primary Component Recreation')
xlabel('Frames')
ylabel('Displacement')
legend('PC1')

%%
clear all;
3

%Loading all data
load('cam1_3.mat')
load('cam2_3.mat')
load('cam3_3.mat')
numFrames1 = size(vidFrames1_3,4);
numFrames2 = size(vidFrames2_3,4);
numFrames3 = size(vidFrames3_3,4);

%%
%Data 1

%Creating Frame for bucket
frame = zeros(480,640);
frame(200:430,300:400) = 1;
frame_u8 = uint8(frame);


%Initializing In-Loop Values
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));


for j = 1:numFrames1
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X1 = vidFrames1_3(:,:,:,j); %Loading Image Frames
    X1_g = rgb2gray(X1); %Converting to Gray
    fX1_g = X1_g .* frame_u8; %Frame only the bucket moving.

    %Focus on the Light by making B/W
    light1 = fX1_g > 250;
    idx1 = find(light1);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light1),idx1);
    
    %Since there are two points, finding the mean point for the can
    mean_x1(j) = mean(X);
    mean_y1(j) = mean(Y);
end

%Inverting image
corrected_image1 = 480 - mean_y1;

%Plotting
figure()
hold on
plot(1:numFrames1,corrected_image1)

%%
%Data 2

%Creating Frame for bucket
frame = zeros(480,640);
frame(100:430,100:460) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));


for j = 1:numFrames2
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X2 = vidFrames2_3(:,:,:,j); %Loading Image Frames
    X2_g = rgb2gray(X2); %Converting to Gray
    fX2_g = X2_g .* frame_u8; %Frame only the bucket moving.

   
    %Focus on the Light by making B/W
    light2 = fX2_g > 240;
    idx2 = find(light2);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light2),idx2);
    
    %Since there are two points, finding the mean point for the can
    mean_x2(j) = mean(X);
    mean_y2(j) = mean(Y);
end

%Inverting image
corrected_image2 = 480 - mean_y2;

%Plotting
plot(1:numFrames2,corrected_image2)

%%
%Data 3

%Creating Frame for bucket
frame = zeros(480,640);
frame(190:330,270:485) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));


for j = 1:numFrames3
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X3 = vidFrames3_3(:,:,:,j); %Loading Image Frames
    X3_g = rgb2gray(X3); %Converting to Gray
    fX3_g = X3_g .* frame_u8; %Frame only the bucket moving.
    
    %Focus on the Light by making B/W
    light3 = fX3_g > 235;
    idx3 = find(light3);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light3),idx3);
    
    %Since there are two points, finding the mean point for the can
    mean_x3(j) = mean(X);
    mean_y3(j) = mean(Y);
end

%Inverting image, using X component as the video is titled
corrected_image3 = 480 - mean_x3;

%Plotting
plot(1:numFrames3,corrected_image3)
title('Uncentered Plots')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')

%%
%All data

%Alligning Data
[M,I1] = min(corrected_image1(1,[1:100]));
min1 = I1;

[M,I2] = min(corrected_image2(1,[1:100]));
min2 = I2;

[M,I3] = min(corrected_image3(1,[1:100]));
min3 = I3;



width = 182;
shift_img1 = 480 - mean_y1(min1:min1+width);
shift_img2 = 480 - mean_y2(min2:min2+width);
shift_img3 = 480 - mean_x3(min3:min3+width);

figure()
hold on
plot(1:width+1,shift_img1)
plot(1:width+1,shift_img2)
plot(1:width+1,shift_img3)
title('Alligned Signals')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')


Full_set = [mean_x1(min1:min1+width); shift_img1; mean_x2(min2:min2+width); shift_img2; shift_img3; mean_y3(min3:min3+width);];
[m,n] = size(Full_set);
mean_FS = mean(Full_set,2);
Full_set = Full_set - repmat(mean_FS,1,n);

% figure()
% plot(1:192,mean_FS)
% title('Mean Signal')
% xlabel('Frames')
% ylabel('Displacement')
% legend('Mean')

[U,S,V] = svd(Full_set/sqrt(n-1));

lambda = diag(S).^2;
sigma = diag(S);
energy = lambda/sum(lambda);


proj = U'*Full_set;

figure()
plot(1:6, [energy(1) sum(energy(1:2)) sum(energy(1:3)) sum(energy(1:4)) sum(energy(1:5)) sum(energy(1:6))],'mo')
title('Energy of each Diagonal Variance')
xlabel('Diagonal Variances')
ylabel('Energy Captured')

figure()
plot(1:width+1,proj(1,:))
title('Primary Component Recreation')
xlabel('Frames')
ylabel('Displacement')
legend('PC1')

%%
clear all;
4

%Loading all data
load('cam1_4.mat')
load('cam2_4.mat')
load('cam3_4.mat')
numFrames1 = size(vidFrames1_4,4);
numFrames2 = size(vidFrames2_4,4);
numFrames3 = size(vidFrames3_4,4);

%%
%Data 1

%Creating Frame for bucket
frame = zeros(480,640);
frame(200:480,300:400) = 1;
frame_u8 = uint8(frame);


%Initializing In-Loop Values
mean_x1 = zeros(1,length(numFrames1));
mean_y1 = zeros(1,length(numFrames1));


for j = 1:numFrames1
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X1 = vidFrames1_4(:,:,:,j); %Loading Image Frames
    X1_g = rgb2gray(X1); %Converting to Gray
    fX1_g = X1_g .* frame_u8; %Frame only the bucket moving.

    %Focus on the Light by making B/W
    light1 = fX1_g > 245;
    idx1 = find(light1);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light1),idx1);
    
    %Since there are two points, finding the mean point for the can
    mean_x1(j) = mean(X);
    mean_y1(j) = mean(Y);
end

%Inverting image
corrected_image1 = 480 - mean_y1;

%Plotting
figure()
hold on
plot(1:numFrames1,corrected_image1)

%%
%Data 2

%Creating Frame for bucket
frame = zeros(480,640);
frame(100:430,100:460) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x2 = zeros(1,length(numFrames2));
mean_y2 = zeros(1,length(numFrames2));


for j = 1:numFrames2
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X2 = vidFrames2_4(:,:,:,j); %Loading Image Frames
    X2_g = rgb2gray(X2); %Converting to Gray
    fX2_g = X2_g .* frame_u8; %Frame only the bucket moving.

   
    %Focus on the Light by making B/W
    light2 = fX2_g > 240;
    idx2 = find(light2);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light2),idx2);
    
    %Since there are two points, finding the mean point for the can
    mean_x2(j) = mean(X);
    mean_y2(j) = mean(Y);
end

%Inverting image
corrected_image2 = 480 - mean_y2;

%Plotting
plot(1:numFrames2,corrected_image2)

%%
%Data 3

%Creating Frame for bucket
frame = zeros(480,640);
frame(190:330,270:485) = 1;
frame_u8 = uint8(frame);

%Initializing In-Loop Values
mean_x3 = zeros(1,length(numFrames3));
mean_y3 = zeros(1,length(numFrames3));


for j = 1:numFrames3
    
    %Process Video so that its in grayscale and framed for just the bucket.
    X3 = vidFrames3_4(:,:,:,j); %Loading Image Frames
    X3_g = rgb2gray(X3); %Converting to Gray
    fX3_g = X3_g .* frame_u8; %Frame only the bucket moving.
    
    %Focus on the Light by making B/W
    light3 = fX3_g > 230;
    idx3 = find(light3);

    %Returning data into Matrix form.
    [Y,X]  = ind2sub(size(light3),idx3);
    
    %Since there are two points, finding the mean point for the can
    mean_x3(j) = mean(X);
    mean_y3(j) = mean(Y);
end

%Inverting image, using X component as the video is titled
corrected_image3 = 480 - mean_x3;

%Plotting
plot(1:numFrames3,corrected_image3)
title('Uncentered Plots')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')

%%
%All data

%Alligning Data
[M,I1] = min(corrected_image1(1,[1:100]));
min1 = I1;

[M,I2] = min(corrected_image2(1,[1:100]));
min2 = I2;

[M,I3] = min(corrected_image3(1,[1:100]));
min3 = I3;



width = 392-I2;
shift_img1 = 480 - mean_y1(min1:min1+width);
shift_img2 = 480 - mean_y2(min2:min2+width);
shift_img3 = 480 - mean_x3(min3:min3+width);

figure()
hold on
plot(1:width+1,shift_img1)
plot(1:width+1,shift_img2)
plot(1:width+1,shift_img3)
title('Alligned Signals')
xlabel('Frames')
ylabel('Displacement')
legend('C1','C2','C3')


Full_set = [mean_x1(min1:min1+width); shift_img1; mean_x2(min2:min2+width); shift_img2; shift_img3; mean_y3(min3:min3+width);];
[m,n] = size(Full_set);
mean_FS = mean(Full_set,2);
Full_set = Full_set - repmat(mean_FS,1,n);

% figure()
% plot(1:192,mean_FS)
% title('Mean Signal')
% xlabel('Frames')
% ylabel('Displacement')
% legend('Mean')

[U,S,V] = svd(Full_set/sqrt(n-1));

lambda = diag(S).^2;
sigma = diag(S);
energy = lambda/sum(lambda);

proj = U'*Full_set;

figure()
plot(1:6, [energy(1) sum(energy(1:2)) sum(energy(1:3)) sum(energy(1:4)) sum(energy(1:5)) sum(energy(1:6))],'mo')
title('Energy of each Diagonal Variance')
xlabel('Diagonal Variances')
ylabel('Energy Captured')

figure()
plot(1:width+1,proj(1,:))
title('Primary Component Recreation')
xlabel('Frames')
ylabel('Displacement')
legend('PC1')
toc
