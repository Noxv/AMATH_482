%Max Walter
%AMATH 482 Test
%Genre Test
clear all; close all; clc
%%
metal1 = load('FFDP.csv');
metal2 = load('MA.csv');
metal3 = load('SK.csv');
metal = [metal1 metal2 metal3];

%%
Pop1 = load('ODE.csv');
Pop2 = load('PM.csv');
Pop3 = load('DP.csv');
pop = [Pop1 Pop2 Pop3];

%%
clas1 = load('JSB.csv');
clas2 = load('LVB.csv');
clas3 = load('WAM.csv');
clas = [clas1 clas2 clas3];
%%
close all;
tot = 1:150;
P = randperm(150,30);
P = sort(P);
% P = 1:5:50;
for i = 1:30
    tot = tot(tot~=P(i));
end
% tot = [1:40];
% P = 41:50;

train_metal = metal(:,tot);
test_metal = metal(:,P);

train_pop = pop(:,tot);
test_pop = pop(:,P);

train_clas = clas(:,tot);
test_clas = clas(:,P);

%%
close all;
feature = 80;
[result,w,U,S,V,threshold1,threshold2] = artist_trianer3(train_metal,train_pop,train_clas,feature);
%[result,w,U,S,V,threshold] = artist_trianer2(train_art1,train_art3,feature);
test_set = [test_metal test_pop test_clas];

TestMat = U'*test_set;
pval = w(5:feature)'*TestMat;
pval = sort(pval);

label = zeros(90,1);
for i = 1:90
    if pval(1,i) < threshold1
        label(i) = 1;
    elseif pval(1,i) > threshold1 && pval(1,i) < threshold2
        label(i) = 3;
    else
        label(i) = 2;
    end
end

hidden = zeros(90,1);
hidden(1:30) = 1;
hidden(31:60) = 2;
hidden(61:90) = 3;
total = 90;

correct = zeros(90,1);
for i = 1:90
    if label(i) == hidden(i)
        correct(i) = 1;
    else
        correct(i) = 0;
    end
end

accuracy = sum(correct)/total;
error = 1 - accuracy;
hold on
plot(pval,1.5,'go')