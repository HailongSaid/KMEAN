function KmeansOfImproved
clear;clc;
global data;
%随机产生数据源
% data1 = [randn(100,2)+ones(100,2);randn(100,2)-ones(100,2)];
load('data1.mat')
k = 6;
figure;
%数据标准化
data = zeros(size(data1));
[data(:,1) me(1) va(1)] = dataNormalization(data1(:,1));
[data(:,2) me(2) va(2)]  = dataNormalization(data1(:,2));

%%
%自己编写的kmeans
[idx c] = kmeansOfMy(data,k);
c = dataRecovery(c,me,va);
%画出各个区域中的散点
count = 0;
for i = 1 : k
    if i == 1
         plot(data1(idx == i,1),data1(idx == 1,2),'r*');
    elseif i == 2
         plot(data1(idx == i,1),data1(idx == i,2),'g*');
    elseif i == 3
         plot(data1(idx == i,1),data1(idx == i,2),'b*');
    elseif i == 4
         plot(data1(idx == i,1),data1(idx == i,2),'ro');
    elseif i == 5
         plot(data1(idx == i,1),data1(idx == i,2),'go');
    elseif i == 6
         plot(data1(idx == i,1),data1(idx == i,2),'bo');
    end
    hold on
    plot(c(i,1),c(i,2),'k^','MarkerSize',16);
    hold on;
    temp = (data(idx == i,1) - c(i,1)) .^ 2 + (data(idx == i,2) - c(i,2)) .^ 2;
    count = count + sum(temp(:));
end
title('My Improved Kmeans');
disp(['My Improved Kmeans中点到各自中心点的距离平方和为：' num2str(count)]);
figure;
%%
%MATLAB自带的函数
[idx,c] = kmeans(data,k);
c = dataRecovery(c,me,va);
%画出各个区域中的散点
count = 0;
for i = 1 : k
    if i == 1
         plot(data1(idx == i,1),data1(idx == 1,2),'r*');
    elseif i == 2
         plot(data1(idx == i,1),data1(idx == i,2),'g*');
    elseif i == 3
         plot(data1(idx == i,1),data1(idx == i,2),'b*');
    elseif i == 4
         plot(data1(idx == i,1),data1(idx == i,2),'ro');
    elseif i == 5
         plot(data1(idx == i,1),data1(idx == i,2),'go');
    elseif i == 6
         plot(data1(idx == i,1),data1(idx == i,2),'bo');
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

%数据复原
function f = dataRecovery(val,me,va)
f = zeros(size(val));
for i = 1 : size(val,2)
    f(:,i) = dataRecovery1(val(:,i),me(i),va(i));
end
end

function f = dataRecovery1(val,me,va)
f = val .* va + me;
end

%数据标准化
function [dataNormal m var] = dataNormalization(val)
m = mean(val);
var = std(val);
dataNormal = (val - m) ./ var;
end

%kmeans函数
function [index c] = kmeansOfMy(data,k)
%kOfVertex = randKOfVertex(k);
kOfVertex = electedInitialCentroid(k);
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
    
%改进处
    %平均距离
    temp = sqrt((vertex(i,1) - data(ix(:),1)) .^ 2 + (vertex(i,2) - data(ix(:),2)) .^ 2);
    s = mean(temp(:));
    
    %找出距离大于2S的点作为该集合的子集
    ixx = [];
    for j = ix
        if sqrt((vertex(i,1) - data(j,1)) .^ 2 + (vertex(i,2) - data(j,2)) .^ 2) > s
            ixx = [ixx j];
        end
    end
    
    %根据该子集重新算出该类的质心
    vertex(i,1:2) = sum(data(ixx(:),:)) ./ length(ixx);
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

%产生k个初始质心
function vertexIndex = electedInitialCentroid(k)
global data;
m = size(data,1);
distans = zeros(m,m);
for i = 1 : m
    distans(i,:) = (data(i,1) - data(:,1)) .^ 2 + (data(i,2) - data(:,2)) .^2;
end
distans = sqrt(distans);
[max1 max2] = find(distans == max(distans(:)));
otherVertexIndex = 1 : m;
vertexIndex1 = [max1(1)  max2(1)];
otherVertexIndex(find(otherVertexIndex - max1(1) == 0)) = [];
otherVertexIndex(find(otherVertexIndex - max2(1) == 0)) = [];
for i = 3 : k
    maxDistans = [0 0];
    for j = otherVertexIndex
        templer = (data(j,1) - data(vertexIndex1,1)) .^ 2 + (data(j,2) - data(vertexIndex1,2)) .^2;
        dis = sum(templer(:));
        if dis >= maxDistans(2)
            maxDistans = [j dis];
        end
    end
    vertexIndex1 = [vertexIndex1 maxDistans(1)];
    otherVertexIndex(find(otherVertexIndex - maxDistans(1) == 0)) = [];
end

vertexIndex = data(vertexIndex1,:);
end
