function KmeansOfOriginal
clear;clc;
global data;
%随机产生数据源
% data = [randn(10,2)+ones(10,2);randn(10,2)-ones(10,2)];
load('data1.mat')
data = data1;
k = 6;
figure;
%%
%自己编写的kmeans
[idx c] = kmeansOfMy(data,k);
%画出各个区域中的散点
count = 0;
for i = 1 : k
    if i == 1
         plot(data(idx == i,1),data(idx == 1,2),'r*');
    elseif i == 2
         plot(data(idx == i,1),data(idx == i,2),'g*');
    elseif i == 3
         plot(data(idx == i,1),data(idx == i,2),'b*');
    elseif i == 4
         plot(data(idx == i,1),data(idx == i,2),'ro');
    elseif i == 5
         plot(data(idx == i,1),data(idx == i,2),'go');
    elseif i == 6
         plot(data(idx == i,1),data(idx == i,2),'bo');
    end
    hold on
    plot(c(i,1),c(i,2),'k^','MarkerSize',16);
    hold on;
    temp = (data(idx == i,1) - c(i,1)) .^ 2 + (data(idx == i,2) - c(i,2)) .^ 2;
    count = count + sum(temp(:));
end
title('My Original Kmeans');
disp('由于这个方法的初始中心是随机选取的，下面的值的每次的取值有可能不一样。');
disp(['My Original Kmeans 中点到各自中心点的距离平方和为：' num2str(count)]);
figure;
%%
%MATLAB自带的函数
[idx,c] = kmeans(data,k);
%画出各个区域中的散点
count = 0;
for i = 1 : k
    if i == 1
         plot(data(idx == i,1),data(idx == 1,2),'r*');
    elseif i == 2
         plot(data(idx == i,1),data(idx == i,2),'g*');
    elseif i == 3
         plot(data(idx == i,1),data(idx == i,2),'b*');
    elseif i == 4
         plot(data(idx == i,1),data(idx == i,2),'ro');
    elseif i == 5
         plot(data(idx == i,1),data(idx == i,2),'go');
    elseif i == 6
         plot(data(idx == i,1),data(idx == i,2),'bo');
    end
    hold on
    plot(c(i,1),c(i,2),'k^','MarkerSize',16);
    hold on;
    temp = (data(idx == i,1) - c(i,1)) .^ 2 + (data(idx == i,2) - c(i,2)) .^ 2;
    count = count + sum(temp(:));
end
title('Matlab Kmeans');
disp(['Matlab Kmeans中点到各自中心点的距离平方和为：' num2str(count)]);
end

%kmeans函数
function [index c] = kmeansOfMy(data,k)
kOfVertex = randElectedInitialCentroid(k);
for i = 1 : size(data,1)
        index(i) = minOfDistans(i,kOfVertex);
end
kOfVertexNew = findVertex(index,k);

while kOfVertexNew ~= kOfVertex
    kOfVertex = kOfVertexNew;
    %将顶点分配到各个区域中去
    for i = 1 : size(data,1)
        index(i) = minOfDistans(i,kOfVertex);
    end

    %找到k个区域的质心
    kOfVertexNew = findVertex(index,k);
end

c = kOfVertex;
end

%找到k个区域的质心
function  vertex = findVertex(index,k)
global data;
vertex = zeros(k,2);
for i = 1 : k
    ix = find(index == i);
    vertex(i,1:2) = sum(data(ix(:),:)) ./ length(ix);
end
end

%找出一个顶点到这些质心的最短的那个质心
function minOfDistans = minOfDistans(vertexIndex,kOfVertex)
global data;
distan = zeros(1,size(kOfVertex,1));
for i = 1 : length(kOfVertex)
    distan(i) = sqrt((data(vertexIndex,1) -  kOfVertex(i,1))^2 + ...
        (data(vertexIndex,2) -  kOfVertex(i,2))^2);
end
minOfDistans = find(distan == min(distan));

end

%随机产生k个质心
function vertexIndex = randElectedInitialCentroid(k)
global data;
vertexIndex1 = zeros(1,k);
vertexIndex1(1) = floor(rand() * size(data,1)) + 1;
count = 1;
while count < k
    num = floor(rand() * size(data,1)) + 1;
   if sum(vertexIndex1 == num) == 0
       vertexIndex1(count + 1) = num;
       count = count + 1;
   end
end

vertexIndex = data(vertexIndex1,:);
end