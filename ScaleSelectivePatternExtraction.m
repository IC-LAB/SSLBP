function SSLBP = ScaleSelectivePatternExtraction(ImageList,Rs,P,ScaleNum,GauFilter,PatternNum,PatternMapping)
% Author: Zhenhua Guo, Xingzheng Wang, Jie Zhou and Jane You
% Paper: Robust texture image representation by scale selective local
% binary patterns, IEEE Transactions on Image Processing
% Version 1.0
% Date: Dec. 16, 2015

% ImageList: image list
% Rs: radius of LBP
% P: neighborhood of LBP
% SclaeNum: parameter to build scale space
% GauFilter: Gaussian smooth filter to build scale space
% 
% PatternMapping: pattern mapping to generate SSLBP

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

for j=1:length(ImageList);
    Gray = imread(ImageList{j});
    Gray = im2double(Gray);
    Gray = (Gray-mean(Gray(:)))/std(Gray(:));
        
    for m=1:length(Rs)
        [CLBP_SC,CLBP_MC] = clbp_single(Gray,Rs(m),P,PatternMapping{m});
        LBPSCTemp{m}(1,:) = CLBP_SC(1:PatternNum)/sum(CLBP_SC);
        LBPMCTemp{m}(1,:) =  CLBP_MC(1:PatternNum)/sum(CLBP_MC);
    end
    for k=2:ScaleNum        
        Gray = imfilter(Gray,GauFilter,'replicate','same');
        for m=1:length(Rs)
            [CLBP_SC,CLBP_MC] = clbp_single(Gray,Rs(m),P,PatternMapping{m});
            LBPSCTemp{m}(k,:) = CLBP_SC(1:PatternNum)/sum(CLBP_SC);
            LBPMCTemp{m}(k,:) =  CLBP_MC(1:PatternNum)/sum(CLBP_MC);
        end
    end
    
    for k=1:length(Rs)
        SSLBP(j,(k-1)*PatternNum+1:k*PatternNum) = max(LBPSCTemp{k});
    end
    
    for k=1:length(Rs)
        SSLBP(j,(k+length(Rs)-1)*PatternNum+1:(k+length(Rs))*PatternNum) = max(LBPMCTemp{k});
    end
    save('tempdata','SSLBP')
end

