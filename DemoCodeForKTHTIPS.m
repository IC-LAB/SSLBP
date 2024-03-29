function DemoCodeForKTHTIPS
% Author: Zhenhua Guo, Xingzheng Wang, Jie Zhou and Jane You
% Paper: Robust texture image representation by scale selective local
% binary patterns, IEEE Transactions on Image Processing
% Version 1.0
% Date: Dec. 16, 2015

rootpic = 'D:\Work Pictures\Texture Database\KTH-TIPS\KTH_TIPS\'; % picture root for the database

% retrieve all image folders of KTH-TIPS database
filelist = dir(rootpic);
foldernum = 0;
for i=3:length(filelist)
    foldername = sprintf('%s%s',rootpic,filelist(i).name);
    if isdir(foldername)
        foldernum = foldernum+1;
        folderlist{foldernum} = foldername;
    end
end

% default parameters
Rs =[3,9];
P = 24;
ScaleNum = 4;
GauFilter = fspecial('gaussian',9,2^(1/4));
PatternNum = 600;

classNum = 10;
picNumPerClass = 81;

% select first 25% images as training set, store filename in ImageList
picCount = 0;
for i=1:classNum
    filelist = dir([folderlist{i},'\*.png']);
    for j=1:round(length(filelist)/4);
        filename = sprintf('%s\\%s', folderlist{i}, filelist(j).name);
        picCount = picCount+1;
        ImageList{picCount} = filename;
    end
end

ScaleSelectivePatterns = GenerateScaleSelectivePatterns(ImageList,Rs,P,ScaleNum,GauFilter,PatternNum);
% save('temp','ScaleSelectivePatterns')

% To speed up feature extraction, create a new mapping for Scale Selective Patterns
% patternMappingri=getmapping(P,'ri');
% save('RIMap24','patternMappingri')
load RIMap24 % it is generated by getmapping(24,'ri')

tabletemp = [patternMappingri.table,patternMappingri.table+patternMappingri.num];
clear patternMappingri

for i=1:length(Rs)        
    patternMappingriNew{i}.num = PatternNum+1;    
    patternMappingriNew{i}.table1 = single(PatternNum*ones(size(tabletemp)));
    patternMappingriNew{i}.table2 = single(PatternNum*ones(size(tabletemp)));
    for k=1:PatternNum
        patternid = ScaleSelectivePatterns(i,k)-1;
        idx = find(tabletemp==patternid);
        patternMappingriNew{i}.table1(idx) = k-1;
        patternid = ScaleSelectivePatterns(i+length(Rs),k)-1;
        idx = find(tabletemp==patternid);
        patternMappingriNew{i}.table2(idx) = k-1;
    end
    patternMappingriNew{i}.samples = P;
end
save('PatternMappingNew','patternMappingriNew')

% "patternMappingriNew" stores new pattern mapping for feature extractions
load PatternMappingNew
% select images as training set, store filename in ImageList
picCount = 0;
for i=1:classNum
    filelist = dir([folderlist{i},'\*.png']);
    for j=1:length(filelist)
        filename = sprintf('%s\\%s', folderlist{i}, filelist(j).name);
        picCount = picCount+1;
        ImageList{picCount} = filename;
    end
end

% Generate scale selective patterns
SSLBP = ScaleSelectivePatternExtraction(ImageList,Rs,P,ScaleNum,GauFilter,PatternNum,patternMappingriNew);    
save('SSLBPFeature','SSLBP')
load SSLBPFeature

trainNumPerClass = 40;
testNum = 100;
% Test on nearest subspace classifier (NSC)
SSLBPHist = SSLBP.^0.5;
Rec = NSC_Random(classNum,picNumPerClass,SSLBPHist,trainNumPerClass,testNum);       

% Test on nearest neighbor classifier (NNC)
for i=1:classNum*picNumPerClass
    FullDist(i,:) = distMATChiSquare(SSLBPHist,SSLBPHist(i,:));
end
Rec=NNC_Random(classNum,picNumPerClass,FullDist,trainNumPerClass,testNum);