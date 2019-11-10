function [HSI_ssd,Sm_HSI_EAS]...
        = Saliency_EAS(scene,varargin)
%extract saliency map from hyperspectral data
%input, scene, hyperspectral image in mat format
%output, 
% HSI_ssd, conspicuity map built from spectral spatial distribution
% Sm_HSI_EAS, saliency map 
% dependence, matlabHyperspectralToolbox
%             L1/2 umixming

clear global
[path, ~] = fileparts(mfilename('fullpath'));

[~,~,b] = size(scene);  
data = normalise(scene, '', 1);
bands = cell(1,b);
for i = 1:b 
    bands{i}=data(:,:,i);
end
% ------------------------------------------------------------------------
% Four spectral groups
% ------------------------------------------------------------------------
bPyr = GaussianPyramid(bands);
% HSI whole spectral bands with Euclidian distance and spectral angle
% distance
% ------------------------------------------------------------------------
% Euclidear distance
FM1 = BandsFeatureMap_2(bPyr);
% Spectral Angle distance
FM2 = BandsFeatureMap_3(bPyr);
HSI_spectralED = ConspicuityMap_2(FM1);
HSI_spectralSAD = ConspicuityMap_2(FM2);
% ------------------------------------------------------------------------
% conspicuity map built from spectral spatial distribution
% ------------------------------------------------------------------------
HSI_ssd =  speSpaDis(data);
% ------------------------------------------------------------------------
% saliency maps generated from different combinations of conspicuity maps
Sm_HSI_EAS = SaliencyMap(HSI_spectralED, HSI_spectralSAD, HSI_ssd);

% ------------------------------------------------------------------------
% Normalize to uint8 0-255
% ------------------------------------------------------------------------

HSI_spectralED=uint8(Normalize(HSI_spectralED)*255);
HSI_spectralSAD=uint8(Normalize(HSI_spectralSAD)*255);
HSI_ssd = uint8(Normalize(HSI_ssd)*255);

Sm_HSI_EAS = uint8(Normalize(Sm_HSI_EAS)*255);

% ------------------------------------------------------------------------
% GaussianPyramid
% ------------------------------------------------------------------------
function pyramid = GaussianPyramid(image)
b = size(image,2);

if iscell(image) == 0
    pyramid=cell(1,9);
    pyramid{1} = image;
    for level=2:9
        im = gausmooth(pyramid{level-1});
        s=ceil(size(pyramid{level-1})/2.0);
        pyramid{level} = imresize(im,s);
    end
else  
    pyramid=cell(b,9);
    for i =1:b
        pyramid{i,1}= image{i};
        for level=2:9;
            im = gausmooth(pyramid{i,level-1});
            s=ceil(size(pyramid{i,level-1})/2.0);
            pyramid{i,level} = imresize(im,s);
        end
    end
end    
% ------------------------------------------------------------------------
% GaussianSmooth
% ------------------------------------------------------------------------
function im = gausmooth(im)
[m,n] = size(im);
GaussianDieOff = .0001;  
pw = 1:30; 
ssq = 2;
width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
if isempty(width)
    width = 1;  % the user entered a really small sigma
end
t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/sum(exp(-(t.*t)/(2*ssq)));     

im = imfilter(im, gau,'conv','replicate');   % run the filter accross rows
im = imfilter(im, gau','conv','replicate'); % and then accross columns


function [rgFM, byFM] = BandsFeatureMap_1(bPyr)
b = size(bPyr,1); % number of bands
g = floor(b/4); % number of bands per group 
for c=3:5
    for delta = 3:4
        s = c + delta;
        rgFM{c,s} = zeros(size(bPyr{1,c}));
        byFM{c,s} = zeros(size(bPyr{1,c}));
        for i = 0:g-1
            rgFM_group = abs(Subtract(bPyr{b-i,c}-bPyr{floor(b/2)-i,c}, bPyr{b-i,s}-bPyr{floor(b/2)-i,s}));
            rgFM{c,s} = rgFM{c,s} + rgFM_group.^2;
            byFM_group = abs(Subtract(bPyr{floor(b*3/4)-i,c}-bPyr{floor(b/4)-i,c}, bPyr{floor(b*3/4)-i,s}-bPyr{floor(b/4)-i,s}));
            byFM{c,s} = byFM{c,s} + byFM_group.^2;      
        end
        rgFM{c,s} = sqrt(rgFM{c,s});
        byFM{c,s} = sqrt(byFM{c,s});
    end
end

function FM = BandsFeatureMap_2(bPyr)

for c=3:5
    for delta = 3:4
        s = c + delta;
        FM{c,s} = abs(Subtract_cell(bPyr(:,c), bPyr(:,s)));
    end
end


function FM = BandsFeatureMap_3(bPyr)

for c=3:5
    for delta = 3:4
        s = c + delta;
        b = size(bPyr,1);
        [mc, nc] = size(bPyr{1,c});
        [ms, ns] = size(bPyr{1,s});
        
        matc=cell2matrix(bPyr(:,c),mc,nc,b);
        mats=cell2matrix(bPyr(:,s),ms,ns,b);
        mats=imresize(mats,[mc,nc],'bilinear');
        N_matc=sqrt(sum(matc.^2,3));
        N_mats=sqrt(sum(mats.^2,3));
        pro=dot(matc,mats,3)./(N_matc.*N_mats);
        FM{c,s} = acos(pro);
    end
end

% ------------------------------------------------------------------------
% ConspicuityMap
% ------------------------------------------------------------------------
function HSI_group = ConspicuityMap_1(varargin)

rgFM = varargin{1};
byFM = varargin{2};

dim=size(rgFM{3,6});

HSI_group=zeros(dim);


    
    for c=3:5
        for delta = 3:4
            s = c + delta;
    %         weight=s*c;
            weight=1;
            HSI_group = Add(HSI_group,weight*Normalize(rgFM{c,s}));
            HSI_group = Add(HSI_group,weight*Normalize(byFM{c,s}));
        end
    end
  
    
function HSI_spectralED = ConspicuityMap_2(varargin)
FM = varargin{1};
dim=size(FM{3,6});

HSI_spectralED=zeros(dim);
for c=3:5
    for delta = 3:4
        s = c + delta;
%         weight=s*c;
        weight=1;
        HSI_spectralED=Add(HSI_spectralED,weight*Normalize(FM{c,s}));
    end
end

% ------------------------------------------------------------------------
% Normalize
% ------------------------------------------------------------------------

function normalized = Normalize(map)

% Normalize map to range [0..1]
minValue = min(min(map));
map = map-minValue;
maxValue = max(max(map));
if maxValue>0
    map = map/maxValue;
end

% Position of local maxima
lmax = LocalMaxima(map);

% Position of global maximum
gmax = (map==1.0);

% Local maxima excluding global maximum
lmax = lmax .* (gmax==0);

% Average of local maxima excluding global maximum
nmaxima=sum(sum(lmax));
if nmaxima>0
    m = sum(sum(map.*lmax))/nmaxima;
else
    m = 0;
end
normalized = map*(1.0-m)^2;

% ------------------------------------------------------------------------
% LocalMaxima
% ------------------------------------------------------------------------

function maxima = LocalMaxima(A)

nRows=size(A,1);
nCols=size(A,2);
% compare with bottom, top, left, right
maxima =           (A > [A(2:nRows, :);   zeros(1, nCols)]);
maxima = maxima .* (A > [zeros(1, nCols); A(1:nRows-1, :)]);
maxima = maxima .* (A > [zeros(nRows, 1), A(:, 1:nCols-1)]);
maxima = maxima .* (A > [A(:, 2:nCols),   zeros(nRows, 1)]);

% ------------------------------------------------------------------------
% Subtract
% ------------------------------------------------------------------------

function result = Subtract(im1, im2)

im2 = imresize(im2, size(im1), 'bilinear');
result = im1 - im2;

% ------------------------------------------------------------------------
% Add
% ------------------------------------------------------------------------

function result = Add(im1, im2)

im2 = imresize(im2, size(im1), 'bilinear');
result = im1 + im2;

function result = Subtract_cell(im1, im2)
b=size(im1,1);
result=zeros(size(im1{1}));
for i=1:b
    im2{i} = imresize(im2{i}, size(im1{i}), 'bilinear');
    result_temp = im1{i} - im2{i};
    result_temp = result_temp.^2;
    result=result+result_temp;  
end

% ------------------------------------------------------------------------
% SaliencyMap
% ------------------------------------------------------------------------

function sm=SaliencyMap(varargin)

iCM = varargin{1};
oCM = varargin{2};

if length(varargin) == 2
    sm=(Normalize(iCM)+Normalize(oCM))/2;
elseif length(varargin) == 3
    cCM = varargin{3};
    sm = (Normalize(iCM)+Normalize(cCM)+Normalize(oCM))/3;
elseif length(varargin) == 4
    sm = (Normalize(iCM)+Normalize(oCM)+Normalize(varargin{3})+Normalize(varargin{4}))/4;    
end

% ------------------------------------------------------------------------
% ShowImage
% ------------------------------------------------------------------------

function ShowImage(image,fTitle)

% figure(nFigure);
figure;
imshow(image);
axis image;
title(fTitle);

function matrix = cell2matrix(cell,m,n,b)
matrix=zeros(m,n,b);
for i=1:b
    matrix(:,:,i)= cell{i,1};
end



