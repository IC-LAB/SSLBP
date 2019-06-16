function Rec = NSC_Random(ClassNum,PicNumPerClass,SSLBPHists,TrainNumPerClass,TestNum)              
% Author: Zhenhua Guo, Xingzheng Wang, Jie Zhou and Jane You
% Paper: Robust texture image representation by scale selective local
% binary patterns, IEEE Transactions on Image Processing
% Version 1.0
% Date: Dec. 16, 2015

% ClassNum: number of class
% PicNumPerClass: number of picture per class
% TrainNumPerClass: number of gallery picture per class
% TestNum: number of test
% SSLBPHists: histogram of SSLLBP

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
    clear Tt_DAT Tr_DAT trls ttls
    Tt_DAT = SSLBPHists(testIDs,:)';
    Tr_DAT = SSLBPHists(trainIDs,:)';
    trls = trainClassIDs;
    ttls = testClassIDs;       
   
    tr_dat  =  Tr_DAT;
    tt_dat  =  Tt_DAT;
    
    tr_dat  =  tr_dat./( repmat(sqrt(sum(tr_dat.*tr_dat)), [size(tr_dat,1),1]) );
    tt_dat  =  tt_dat./( repmat(sqrt(sum(tt_dat.*tt_dat)), [size(tr_dat,1),1]) );
    
    for j=1:ClassNum
        tridx = find(trls==j);        
        s=pinv(tr_dat(:,tridx))*tt_dat;
        err(j,:)=sqrt(sum((tr_dat(:,tridx)*s-tt_dat).^2,1));
    end
    [distNew, ID]= min(err);    
    cornum      =   sum(ID==ttls);
    Rec(k)         =   [cornum/length(ttls)]*100; % recognition rate
end