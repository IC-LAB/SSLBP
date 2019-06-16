function [CLBP_SC,CLBP_MC] = clbp_single(varargin) % image,radius,neighbors,mapping,mode)
% Author: Zhenhua Guo, Xingzheng Wang, Jie Zhou and Jane You
% Paper: Robust texture image representation by scale selective local
% binary patterns, IEEE Transactions on Image Processing
% Version 1.0
% Date: Dec. 16, 2015
% Acknowledge: part of codes are copied from MVG of Oulu University

% image: image list of training sample
% Rs: radius of LBP
% P: neighborhood of LBP
% mapping: pattern mapping
% mode: generate pattern histogram or pattern index for each pixel

% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};

d_image=double(image);
if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);

    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
end

% Determine the dimensions of the input image.
[ysize xsize] = size(image);

miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C);
bins = 2^neighbors;

% Initialize the result matrix with zeros.
CLBP_SC=zeros(dy+1,dx+1);
CLBP_MC=zeros(dy+1,dx+1);
CLBP_C=zeros(dy+1,dx+1);

% to speed up feature extraction, compute "single" float here
spoints = single(spoints);
d_C = single(C*100);
d_image=single(image*100);
CLBP_SC=single(CLBP_SC);
CLBP_MC=single(CLBP_MC);

%Compute the LBP code image
for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = d_image(ry:ry+dy,rx:rx+dx);
    D{i} = N >= d_C;   
    Diff{i} = abs(N-d_C);    
    MeanDiff(i) = mean(mean(Diff{i}));
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = (1 - tx) * (1 - ty);
    w2 =      tx  * (1 - ty);
    w3 = (1 - tx) *      ty ;
    w4 =      tx  *      ty ;
    w1 = single(w1);
    w2 = single(w2);
    w3 = single(w3);
    w4 = single(w4);
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    D{i} = N >= d_C;  
    
    Diff{i} = abs(N-d_C);    
    MeanDiff(i) = mean(mean(Diff{i}));
  end  
end
% Difference threshold for CLBP_M
DiffThreshold = mean(MeanDiff);
% compute CLBP_S and CLBP_M
for i=1:neighbors
  % Update the result matrix.
  v = 2^(i-1);
  CLBP_SC = CLBP_SC + v*D{i};
  CLBP_MC = CLBP_MC + v*(Diff{i}>=DiffThreshold);
end
% CLBP_C
CLBP_C = d_C>=mean(d_image(:));
CLBP_SC = CLBP_C*2^24+double(CLBP_SC); % 2^24 single will bring an error
CLBP_MC = CLBP_C*2^24+double(CLBP_MC);

sizarray = size(CLBP_SC);
CLBP_SC = CLBP_SC(:);
CLBP_MC = CLBP_MC(:);
CLBP_SC = mapping.table1(CLBP_SC+1);
CLBP_MC = mapping.table2(CLBP_MC+1);

% pixel or histogram
if (strcmp(mode,'h'))
    CLBP_SC = hist(CLBP_SC,mapping.num);
    CLBP_MC = hist(CLBP_MC,mapping.num);
else
    CLBP_SC = reshape(CLBP_SC,sizarray);
    CLBP_MC = reshape(CLBP_MC,sizarray);
end