function [result,w,U,S,V,threshold1,threshold2] = artist_trianer3(t_art1,t_art2,t_art3,feature)

nart1 = length(t_art1(1,:));
nart2 = length(t_art2(1,:));
nart3 = length(t_art3(1,:));


[U,S,V] = svd([t_art1 t_art2 t_art3],'econ');

artists = S*V';
U = U(:, 5:feature);

t_art1 = artists(1:feature,1:nart1);
t_art2 = artists(1:feature,nart1+1:nart1+nart2);
t_art3 = artists(1:feature,(nart1+nart2)+1:nart1+nart2+nart3);
art_all = [t_art1 t_art2 t_art3];

mart1 = mean(t_art1,2);
mart2 = mean(t_art2,2);
mart3 = mean(t_art3,2);
mean_all = mean(art_all,2);

Sw = 0;

for i = 1:nart1
    Sw = Sw + (t_art1(:,i)-mart1)*(t_art1(:,i)-mart1)';
end

for i = 1:nart2
    Sw = Sw + (t_art2(:,i)-mart2)*(t_art2(:,i)-mart2)';
end


for i = 1:nart3
    Sw = Sw + (t_art3(:,i)-mart3)*(t_art3(:,i)-mart3)';
end



Sb = (mart1 - mean_all)*(mart1 - mean_all)' + (mart2 - mean_all)*(mart2 - mean_all)' + (mart3 - mean_all)*(mart3 - mean_all)';

Sb = Sb/3;


[V2,D] = eig(Sb,Sw);
[lambda,ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

vart1 = w'*t_art1;
vart2 = w'*t_art2;
vart3 = w'*t_art3;
result = [vart1,vart2,vart3];

sortart1 = sort(vart1);
sortart2 = sort(vart2);
sortart3 = sort(vart3);

t1 = length(sortart1);
t2 = 1;

while sortart1(t1)>sortart3(t2)
    t1 = t1 - 1;
    t2 = t2+1;

end
threshold1 = (sortart3(t1)+sortart1(t2))/2;

t2 = length(sortart1);
t3 = 1;

while sortart1(t2)>sortart2(t3)
    t2 = t2 - 1;
    t3 = t3+1;
end


threshold2 = (sortart2(t2)+sortart3(t3))/2;

plot(sortart1,0,'bo')
hold on
plot(sortart2,1,'ro')
plot(sortart3,2,'ko')
title('Blue: Artist 1, Red: Artist 2, Black: Artist 3, Green: Test Data')
xlabel('Projected Axis');
plot([threshold1 threshold1],[0 2.5],'r')
plot([threshold2 threshold2],[0 2.5],'b')