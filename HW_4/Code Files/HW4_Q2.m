%Max Walter
%AMATH 482 Test
%Same Genre TEst
clear all; close all; clc
%%
art1 = load('FFDP.csv');
%%
art2 = load('MA.csv');
%%
art3 = load('SK.csv');

%%
close all;
tot = 1:50;
P = randperm(50,20);
P = sort(P);
% P = 1:5:50;
for i = 1:10
    tot = tot(tot~=P(i));
end
% tot = [1:40];
% P = 41:50;

train_art1 = art1(:,tot);
test_art1 = art1(:,P);

train_art2 = art2(:,tot);
test_art2 = art2(:,P);

train_art3 = art3(:,tot);
test_art3 = art3(:,P);

close all;
feature = 40;
[result,w,U,S,V,threshold1,threshold2] = artist_trianer3(train_art1,train_art2,train_art3,feature);
%[result,w,U,S,V,threshold] = artist_trianer2(train_art1,train_art3,feature);
test_set = [test_art1 test_art2 test_art3];

TestMat = U'*test_set;
pval = w(5:feature)'*TestMat;
pval = sort(pval);

label = zeros(30,1);
for i = 1:30
    if pval(1,i) < threshold1
        label(i) = 1;
    elseif pval(1,i) > threshold1 && pval(1,i) < threshold2
        label(i) = 3;
    else
        label(i) = 2;
    end
end

hidden = zeros(30,1);
hidden(1:10) = 1;
hidden(11:20) = 2;
hidden(21:30) = 3;
total = 30;

correct = zeros(30,1);
for i = 1:30
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