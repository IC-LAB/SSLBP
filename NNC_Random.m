function Rec=NNC_Random(ClassNum,PicNumPerClass,FullDist,TrainNumPerClass,TestNum)
% Author: Zhenhua Guo, Xingzheng Wang, Jie Zhou and Jane You
% Paper: Robust texture image representation by scale selective local
% binary patterns, IEEE Transactions on Image Processing
% Version 1.0
% Date: Dec. 16, 2015

% ClassNum: number of class
% PicNumPerClass: number of picture per class
% FullDist: full distance matrix
% TestNum: number of test
% TrainNumPerClass: number of gallery picture per class

if nargin<4
    TrainNumPerClass = 40;
end

if nargin<5
    TestNum = 100;
end
     
for k=1:TestNum
    p = randperm(PicNumPerClass);
    randTrain = p(1:TrainNumPerClass);
    randTest = p(TrainNumPerClass+1:PicNumPerClass);
    for i=1:ClassNum
        trainIDs((i-1)*TrainNumPerClass+1:i*TrainNumPerClass) = randTrain+(i-1)*PicNumPerClass;
        trainClassIDs((i-1)*TrainNumPerClass+1:i*TrainNumPerClass) = i;
        
        testIDs((i-1)*(PicNumPerClass-TrainNumPerClass)+1:i*(PicNumPerClass-TrainNumPerClass)) = randTest+(i-1)*PicNumPerClass;
        testClassIDs((i-1)*(PicNumPerClass-TrainNumPerClass)+1:i*(PicNumPerClass-TrainNumPerClass)) = i;
    end
    
    Dist = FullDist(:,trainIDs);
    Dist = Dist(testIDs,:);
    [distNew, index]= min(Dist,[],2);
    idx = find(trainClassIDs(index)==testClassIDs);
    Rec(k)=length(idx)/length(testIDs)*100; % recognition rate
end