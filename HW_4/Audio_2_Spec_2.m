%AMATH 482-HW 4
%Audio to Spec File 2

data_path = 'C:\Users\maxim\MATLAB_DRIVE\Class\AMATH\AMATH_482\HW_4';
audio_file = '\Audio_Files2';
artist = '\JSB';
prefix = '\JSB_';
ext = '.m4a';
FS = 44100;
spec_out = [];
step = 0;
run = 0;

for j = 1:10
    clear data
    dir = [data_path audio_file artist prefix num2str(j) ext];
    data = audioread(dir);
    %data = downsample(data,10);
    data = data';
    leng = FS * 5;
    step = 0;
    for i = 1:5
        file = data(1,(step + 1):(i*leng));
        spec = spectrogram(file);
        spec_out(:,i + run) = spec(:);
        step = step + leng;
    end
    run = run + 5;
end

writematrix(spec_out,'JSB.csv')