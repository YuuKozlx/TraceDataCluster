clear;clc;close

load matchResult.txt;

figure;
hold on
for i = 1:12
plot(matchResult(:,2*i - 1),matchResult(:,2*i));
end


