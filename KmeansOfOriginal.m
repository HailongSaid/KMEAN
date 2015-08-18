function KmeansOfOriginal
clear;clc;
global data;
%�����������Դ
% data = [randn(10,2)+ones(10,2);randn(10,2)-ones(10,2)];
load('data1.mat')
data = data1;
k = 6;
figure;
%%
%�Լ���д��kmeans
[idx c] = kmeansOfMy(data,k);
%�������������е�ɢ��
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
disp('������������ĳ�ʼ���������ѡȡ�ģ������ֵ��ÿ�ε�ȡֵ�п��ܲ�һ����');
disp(['My Original Kmeans �е㵽�������ĵ�ľ���ƽ����Ϊ��' num2str(count)]);
figure;
%%
%MATLAB�Դ��ĺ���
[idx,c] = kmeans(data,k);
%�������������е�ɢ��
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
disp(['Matlab Kmeans�е㵽�������ĵ�ľ���ƽ����Ϊ��' num2str(count)]);
end

%kmeans����
function [index c] = kmeansOfMy(data,k)
kOfVertex = randElectedInitialCentroid(k);
for i = 1 : size(data,1)
        index(i) = minOfDistans(i,kOfVertex);
end
kOfVertexNew = findVertex(index,k);

while kOfVertexNew ~= kOfVertex
    kOfVertex = kOfVertexNew;
    %��������䵽����������ȥ
    for i = 1 : size(data,1)
        index(i) = minOfDistans(i,kOfVertex);
    end

    %�ҵ�k�����������
    kOfVertexNew = findVertex(index,k);
end

c = kOfVertex;
end

%�ҵ�k�����������
function  vertex = findVertex(index,k)
global data;
vertex = zeros(k,2);
for i = 1 : k
    ix = find(index == i);
    vertex(i,1:2) = sum(data(ix(:),:)) ./ length(ix);
end
end

%�ҳ�һ�����㵽��Щ���ĵ���̵��Ǹ�����
function minOfDistans = minOfDistans(vertexIndex,kOfVertex)
global data;
distan = zeros(1,size(kOfVertex,1));
for i = 1 : length(kOfVertex)
    distan(i) = sqrt((data(vertexIndex,1) -  kOfVertex(i,1))^2 + ...
        (data(vertexIndex,2) -  kOfVertex(i,2))^2);
end
minOfDistans = find(distan == min(distan));

end

%�������k������
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