function ScaleSelectivePatterns = GenerateScaleSelectivePatterns(ImageList,Rs,P,ScaleNum,GauFilter,PatternNum)
% Author: Zhenhua Guo, Xingzheng Wang, Jie Zhou and Jane You
% Paper: Robust texture image representation by scale selective local
% binary patterns, IEEE Transactions on Image Processing
% Version 1.0
% Date: Dec. 16, 2015

% ImageList: image list of training sample
% Rs: radius of LBP
% P: neighborhood of LBP
% SclaeNum: parameter to build scale space
% GauFilter: Gaussian smooth filter to build scale space
% PatternNum: number of scale selective patterns

% default parameters
if nargin<2
    Rs = [3,9];
end
if nargin<3
    P = 24;
end
if nargin<4
    ScaleNum = 4;
end
if nargin<5
    GauFilter = fspecial('gaussian',9,2^(1/4));
end
if nargin<6
    PatternNum = 600;
end

% generate mapping
patternMappingri=getmapping(P,'ri');
% load RIMap24 % it is generated by getmapping(24,'ri')

LBPSigCenterSum = zeros(length(Rs),patternMappingri.num*2);
LBPMagCenterSum = zeros(length(Rs),patternMappingri.num*2);

for i=1:length(ImageList)
    Gray = imread(ImageList{i});
    Gray = im2double(Gray);
    Gray = (Gray-mean(Gray(:)))/std(Gray(:));    % get 0-mean and 1-standard
    
    for m=1:length(Rs)
        [CLBP_S,CLBP_M,CLBP_C] = clbp(Gray,Rs(m),P,patternMappingri,'x');        
        CLBP_SC = ([CLBP_S(:),CLBP_C(:)]);
        Hist3D = hist3(CLBP_SC,[patternMappingri.num,2]);
        LBPSigCenterTemp{m}(1,:) = reshape(Hist3D,1,numel(Hist3D))/numel(CLBP_S);
        CLBP_MC = ([CLBP_M(:),CLBP_C(:)]);
        Hist3D = hist3(CLBP_MC,[patternMappingri.num,2]);
        LBPMagCenterTemp{m}(1,:) = reshape(Hist3D,1,numel(Hist3D))/numel(CLBP_M);
    end
    for k=2:ScaleNum        
        Gray = imfilter(Gray,GauFilter,'replicate','same');
        for m=1:length(Rs)
            [CLBP_S,CLBP_M,CLBP_C] = clbp(Gray,Rs(m),P,patternMappingri,'x');            
            CLBP_SC = ([CLBP_S(:),CLBP_C(:)]);
            Hist3D = hist3(CLBP_SC,[patternMappingri.num,2]);
            LBPSigCenterTemp{m}(k,:) = reshape(Hist3D,1,numel(Hist3D))/numel(CLBP_S);
            CLBP_MC = ([CLBP_M(:),CLBP_C(:)]);
            Hist3D = hist3(CLBP_MC,[patternMappingri.num,2]);
            LBPMagCenterTemp{m}(k,:) = reshape(Hist3D,1,numel(Hist3D))/numel(CLBP_M);
        end
    end
    % scale selective
    for m=1:length(Rs)
        LBPSigCenterSum(m,:) = LBPSigCenterSum(m,:)+max(LBPSigCenterTemp{m});
        LBPMagCenterSum(m,:) = LBPMagCenterSum(m,:)+max(LBPMagCenterTemp{m});
    end
end

% get pattern index
for m=1:length(Rs)    
    [LBPTemp,patternidxall] = sort(LBPSigCenterSum(m,:),'descend');
    ScaleSelectivePatterns(m,1:PatternNum) = patternidxall(1:PatternNum);
end
for m=1:length(Rs)    
    [LBPTemp,patternidxall] = sort(LBPMagCenterSum(m,:),'descend');
    ScaleSelectivePatterns(m+length(Rs),1:PatternNum) = patternidxall(1:PatternNum);
end