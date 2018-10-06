function [loc, lev] = SubFun_judgeLev(X,judgePoint)
% X,矢量原始数据；
% judgePoint 判断间断点，如6个level需要6个判断点；
lenT = length(X);
lenJ = length(judgePoint);
lev = zeros(lenT,1);
loc = zeros(lenT,1);
LocIntervel = judgePoint-circshift([judgePoint(1:end-1) 0],1);
for i = 1:lenT
    for j = 1:lenJ
        if X(i)<=judgePoint(j)
           lev(i) =  j;
             if j==1
             loc(i) = X(i)/LocIntervel(1);
             else
              loc(i) = (X(i)-judgePoint(j-1))/LocIntervel(j-1);   
             end
           break;
        end
    end
end