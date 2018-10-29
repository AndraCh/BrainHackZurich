%% Andra Chincisan, Institute of Neuropathology, USZ
% 2018
%%
close all;
clc;
clear all;
%% Read image
Cmatrix=[];
image_type = '.tif';
path_dir = 'Images\'
path_dir_masks = 'Images_mask\'
srcFiles = dir(strcat(path_dir, '*', image_type));
for i = 1: length(srcFiles)
    image_no = i;
    close all;
    % Read image
    image = imread (strcat(path_dir,srcFiles(i).name)); 
    % Gray image
    imageGray = getgrayimage (image);
    % Segment image
    imageBinary = segmentimage (imageGray);
    figure, imshow(image(:,:,3)); title('imageVacoulesNew');
    % Get vacuoles
    imageVacoules = getvacoulesimage (imageBinary);
    % Process vacuoles
    imageVacoulesNew= remove_vessels (imageGray, imageVacoules);
    % Quantify vacuoles
    [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules, Area_all,Eccentricity_all] = quantifyvacoules (imageVacoulesNew);
    figure, imshow(imageVacoulesNew); title('imageVacoulesNew');
    % Fuse original and segmented image
    C = imfuse(imageGray,imageVacoulesNew);
    figure, imshow(C); title(strcat(path_dir,srcFiles(i).name));
    % Save fused image
    %imwrite(mat2gray(C),strcat(path_save, srcFiles(i).name));
    % Create arry with data for each each
    %name_image =string (srcFiles(i).name);
   % all_param(i, :) = [{char(name_image)} numberOfVacoules sumAreaVacoules sumPerimeterVacoules Area_all]
    [matrix_features] = compute_features (imageVacoulesNew, imageGray, path_dir_masks);
    labels = check_labels (path_dir_masks, image_no, imageVacoulesNew);
    matrix_features = [matrix_features labels'];
    final_mat = matrix_features;
   
    Cmatrix = vertcat(Cmatrix,final_mat);
end

xlswrite('features.xlsx',Cmatrix) % Save data in Excel file

function labels = check_labels (path_dir_masks, image_no, imageBW)
srcFiles = dir(strcat(path_dir_masks, '*','.tif'));
img = imread (strcat(path_dir_masks,srcFiles(image_no).name)); 
disp(strcat(path_dir_masks,srcFiles(image_no).name))

mask = im2bw(img,0.2);
mask = imcomplement(mask);
mask = bwareaopen(mask,10);
figure, imshow (mask);
[imageLabelM, numberOfObjMask] = bwlabel(mask);
[imageLabelBW, numberOfObjBW] = bwlabel(imageBW);
labels = zeros(1,numberOfObjBW);
for i=1:numberOfObjBW
    for j=1:numberOfObjMask
        if nnz((imageLabelBW==i) & (imageLabelM==j)) > 0
            labels(i) = 1;
        end
    end
end
end

%%
function [matrix_features] = compute_features (image, imageGray, path_dir_masks)
%[left, top, width, height]
imageCC = bwconncomp(image, 4);
S = regionprops(image, 'BoundingBox'); 
features =[S.BoundingBox];
[size_m, size_n] = size(image);
all_the_x1s = features(1:4:end);
all_the_y1s = features(2:4:end);
all_the_widths =  features(3:4:end);
all_the_heights = features(4:4:end);

% Determine the x2's and y2's.
all_the_x2s = all_the_x1s + all_the_widths;
all_the_y2s = all_the_y1s + all_the_heights;

matrix_features = zeros(length(S),14);
for p=1:length(S)
    img = zeros(all_the_heights(p),all_the_widths(p));%all_the_heights(1));
    x=1;
    for i=round(all_the_y1s(p)): round(all_the_y2s(p)) %round(all_the_x1s(1) + all_the_widths(1))
        y = 1;
        for j=round( all_the_x1s(p)) : round(all_the_x2s (p)) %round(all_the_y1s(1) +all_the_heights(1) )
            if i > size_m i = size_m; 
            end
            if j > size_n j = size_n ;
            end
            img(x,y) = imageGray(i,j);
            y = y + 1;
        end
        x = x + 1;
    end
    img = mat2gray(img);
    glcm = graycomatrix(img, 'offset', [0 1], 'Symmetric', true);
    x = haralickTextureFeatures(glcm);
    matrix_features(p,:) = x;
    figure, imshow (img); title ('abc');
    clear img;
end

%labels = check_labels (path_dir_masks, image_no, image);

end
%% Get second layer of the RBG image (blue layer) as gray image
function imageGray = getgrayimage (image)
imageGray = image(:,:,2);
%figure, imshow (imageGray); title ('Gray image'); %imcontrast
end 

%% Threshold segmentation
function imageBinary = segmentimage (imageGray)
% Adjust image contrast
%imageGray = imadjust(imageGray,[0.2 0.8],[]);  
% Binarize image
imageBinary = im2bw(imageGray,0.74);
figure, imshow (imageBinary); title('imageBinary');
% Remove small objects
imageBinary = bwareaopen(imageBinary,20);
%Dilation
se = strel('disk',2);
imageBinary = imdilate(imageBinary,se);
% Erosion
se = strel('disk',2);
imageBinary = imdilate(imageBinary,se);
imageBinary = imfill(imageBinary,'holes');
%figure, imshow(imageBinary);
end

%% Get vacoules
function imageVacoules = getvacoulesimage (imageBinary)
%  Connected components
imageCC = bwconncomp(imageBinary, 4 );
% Geometrical parameters of the image
S = regionprops(imageCC, 'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter'); 
% Label objects in the image
imageLabel = labelmatrix(imageCC);
% Calculate object's shape
shape = 4*pi*[S.Area]./power ([S.Perimeter],2);
% Select smaller and rounder objects
imgBinaryObjSmall = ismember(imageLabel, find(([S.Area] >= 200) & ([S.Area]< 2000) &(shape>= 0.3)));% & ([S.Eccentricity]>0.5)));
% Combine images
imageVacoules = imgBinaryObjSmall;
%figure, imshow (imageVacoules); title ('Final vacoules images');
end

%% Final objects list
function [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules,mean_area_all,mean_eccentricity_all] = quantifyvacoules (imageVacoules)
[imageLabel, numberOfVacoules] = bwlabel(imageVacoules);
 imageCC = bwconncomp(imageVacoules, 4 );
% Geometrical parameters of the image
S = regionprops(imageCC, 'Area', 'Perimeter', 'Eccentricity'); 
% Calculate sum of all vacuole areas in an image
sumAreaVacoules = sum([S.Area]);
% Calculate sum of all vacuole perimeters in an image
sumPerimeterVacoules = sum([S.Perimeter]);
% Calculate mean area of vacuole in an image
mean_area_all = mean([S.Area]);
% Calculate mean eccentricity of vacuole in an image
mean_eccentricity_all= mean([S.Eccentricity]);
end

function mat = window_image(Cx,Cy,D,imageGray)
[max_n, max_m] = size(imageGray);
D = D/3;
new_coordinates = [round(Cx-D/2),round(Cx+D/2),round(Cy-D/2),round(Cy+D/2)];
for i=1:length(new_coordinates)
    if new_coordinates(i) <=0
        new_coordinates(i) = 0;  
    end
end
if (Cx-D)>=0 & (Cx-D)<=max_m
    new_coordinates(1) = round(Cx-D);
end
if (Cx+D)>=0 & (Cx+D)<=max_m
    new_coordinates(2) = round(Cx+D);
end
if (Cy-D)>=0 & (Cy-D)<=max_n
    new_coordinates(3) = round(Cy-D);
end
if (Cy+D)>=0 & (Cy+D)<=max_n
    new_coordinates(4) = round(Cy+D);
end
mat = zeros ( new_coordinates(4)-new_coordinates(3),new_coordinates(2)-new_coordinates(1));
for i=1:(new_coordinates(4)-new_coordinates(3))
    for j=1:(new_coordinates(2)-new_coordinates(1))
        if (new_coordinates(3)+i<= max_n) & (new_coordinates(1)+j<= max_m)
            mat(i,j) = imageGray(new_coordinates(3)+i,new_coordinates(1)+j);
        end
    end
end
end

function [imageVacoules]= remove_vessels (imageGray, imageVacoules)
% Create a small image for each of the initial objects that are considered
% vacuoles. If dark regions appear in the close proximity of an a object it 
% means that there is a vessel close to the object and that this object is 
% not a real vacuole, so it should be remove from the vacuoles list.
% Compare the average intensity of small images with a threshold T =
% 150. Objects with an average intensity <150 will be removed.
 
% Image parameters
S = regionprops(imageVacoules, 'centroid', 'MajorAxisLength'); 
centroids = cat(1, S.Centroid);
diameters = [S.MajorAxisLength];
[imgLabel, numberOfObject] = bwlabel(imageVacoules);

[s1, s2] = size(imageVacoules);
for noObj=1:numberOfObject
     mat = window_image(centroids(noObj,1),centroids(noObj,2),diameters(noObj),imageGray);
     [m, n] = size(mat);
     if (m>10) & (n>10)
        if mean(mean(mat))<150
            imageVacoules(imgLabel==noObj) = 0;
        end
     end
end
% Count objects
[imgLabel, numberOfObject] = bwlabel(imageVacoules);
end

function [out] = GLCM_Features1(glcmin,pairs)
% 
% GLCM_Features1 helps to calculate the features from the different GLCMs
% that are input to the function. The GLCMs are stored in a i x j x n
% matrix, where n is the number of GLCMs calculated usually due to the
% different orientation and displacements used in the algorithm. Usually
% the values i and j are equal to 'NumLevels' parameter of the GLCM
% computing function graycomatrix(). Note that matlab quantization values
% belong to the set {1,..., NumLevels} and not from {0,...,(NumLevels-1)}
% as provided in some references
% http://www.mathworks.com/access/helpdesk/help/toolbox/images/graycomatrix
% .html
% 
% Although there is a function graycoprops() in Matlab Image Processing
% Toolbox that computes four parameters Contrast, Correlation, Energy,
% and Homogeneity. The paper by Haralick suggests a few more parameters
% that are also computed here. The code is not fully vectorized and hence
% is not an efficient implementation but it is easy to add new features
% based on the GLCM using this code. Takes care of 3 dimensional glcms
% (multiple glcms in a single 3D array)
% 
% If you find that the values obtained are different from what you expect 
% or if you think there is a different formula that needs to be used 
% from the ones used in this code please let me know. 
% A few questions which I have are listed in the link 
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/239608
%
% I plan to submit a vectorized version of the code later and provide 
% updates based on replies to the above link and this initial code. 
%
% Features computed 
% Autocorrelation: [2]                      (out.autoc)
% Contrast: matlab/[1,2]                    (out.contr)
% Correlation: matlab                       (out.corrm)
% Correlation: [1,2]                        (out.corrp)
% Cluster Prominence: [2]                   (out.cprom)
% Cluster Shade: [2]                        (out.cshad)
% Dissimilarity: [2]                        (out.dissi)
% Energy: matlab / [1,2]                    (out.energ)
% Entropy: [2]                              (out.entro)
% Homogeneity: matlab                       (out.homom)
% Homogeneity: [2]                          (out.homop)
% Maximum probability: [2]                  (out.maxpr)
% Sum of sqaures: Variance [1]              (out.sosvh)
% Sum average [1]                           (out.savgh)
% Sum variance [1]                          (out.svarh)
% Sum entropy [1]                           (out.senth)
% Difference variance [1]                   (out.dvarh)
% Difference entropy [1]                    (out.denth)
% Information measure of correlation1 [1]   (out.inf1h)
% Informaiton measure of correlation2 [1]   (out.inf2h)
% Inverse difference (INV) is homom [3]     (out.homom)
% Inverse difference normalized (INN) [3]   (out.indnc) 
% Inverse difference moment normalized [3]  (out.idmnc)
%
% The maximal correlation coefficient was not calculated due to
% computational instability 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
%
% Formulae from MATLAB site (some look different from
% the paper by Haralick but are equivalent and give same results)
% Example formulae: 
% Contrast = sum_i(sum_j(  (i-j)^2 * p(i,j) ) ) (same in matlab/paper)
% Correlation = sum_i( sum_j( (i - u_i)(j - u_j)p(i,j)/(s_i.s_j) ) ) (m)
% Correlation = sum_i( sum_j( ((ij)p(i,j) - u_x.u_y) / (s_x.s_y) ) ) (p[2])
% Energy = sum_i( sum_j( p(i,j)^2 ) )           (same in matlab/paper)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + |i-j|) ) ) (as in matlab)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + (i-j)^2) ) ) (as in paper)
% 
% Where:
% u_i = u_x = sum_i( sum_j( i.p(i,j) ) ) (in paper [2])
% u_j = u_y = sum_i( sum_j( j.p(i,j) ) ) (in paper [2])
% s_i = s_x = sum_i( sum_j( (i - u_x)^2.p(i,j) ) ) (in paper [2])
% s_j = s_y = sum_i( sum_j( (j - u_y)^2.p(i,j) ) ) (in paper [2])
%
% 
% Normalize the glcm:
% Compute the sum of all the values in each glcm in the array and divide 
% each element by it sum
%
% Haralick uses 'Symmetric' = true in computing the glcm
% There is no Symmetric flag in the Matlab version I use hence
% I add the diagonally opposite pairs to obtain the Haralick glcm
% Here it is assumed that the diagonally opposite orientations are paired
% one after the other in the matrix
% If the above assumption is true with respect to the input glcm then
% setting the flag 'pairs' to 1 will compute the final glcms that would result 
% by setting 'Symmetric' to true. If your glcm is computed using the
% Matlab version with 'Symmetric' flag you can set the flag 'pairs' to 0
%
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973
% 2. L. Soh and C. Tsatsoulis, Texture Analysis of SAR Sea Ice Imagery
% Using Gray Level Co-Occurrence Matrices, IEEE Transactions on Geoscience
% and Remote Sensing, vol. 37, no. 2, March 1999.
% 3. D A. Clausi, An analysis of co-occurrence texture statistics as a
% function of grey level quantization, Can. J. Remote Sensing, vol. 28, no.
% 1, pp. 45-62, 2002
% 4. http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
%
%
% Example:
%
% Usage is similar to graycoprops() but needs extra parameter 'pairs' apart
% from the GLCM as input
% I = imread('circuit.tif');
% GLCM2 = graycomatrix(I,'Offset',[2 0;0 2]);
% stats = GLCM_features1(GLCM2,0)
% The output is a structure containing all the parameters for the different
% GLCMs
%
% [Avinash Uppuluri: avinash_uv@yahoo.com: Last modified: 11/20/08]
% If 'pairs' not entered: set pairs to 0 
if ((nargin > 2) || (nargin == 0))
   error('Too many or too few input arguments. Enter GLCM and pairs.');
elseif ( (nargin == 2) ) 
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
elseif (nargin == 1) % only GLCM is entered
    pairs = 0; % default is numbers and input 1 for percentage
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
       error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
       error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end    
end
format long e
if (pairs == 1)
    newn = 1;
    for nglcm = 1:2:size(glcmin,3)
        glcm(:,:,newn)  = glcmin(:,:,nglcm) + glcmin(:,:,nglcm+1);
        newn = newn + 1;
    end
elseif (pairs == 0)
    glcm = glcmin;
end
size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);
% checked 
out.autoc = zeros(1,size_glcm_3); % Autocorrelation: [2] 
out.contr = zeros(1,size_glcm_3); % Contrast: matlab/[1,2]
out.corrm = zeros(1,size_glcm_3); % Correlation: matlab
out.corrp = zeros(1,size_glcm_3); % Correlation: [1,2]
out.cprom = zeros(1,size_glcm_3); % Cluster Prominence: [2]
out.cshad = zeros(1,size_glcm_3); % Cluster Shade: [2]
out.dissi = zeros(1,size_glcm_3); % Dissimilarity: [2]
out.energ = zeros(1,size_glcm_3); % Energy: matlab / [1,2]
out.entro = zeros(1,size_glcm_3); % Entropy: [2]
out.homom = zeros(1,size_glcm_3); % Homogeneity: matlab
out.homop = zeros(1,size_glcm_3); % Homogeneity: [2]
out.maxpr = zeros(1,size_glcm_3); % Maximum probability: [2]
out.sosvh = zeros(1,size_glcm_3); % Sum of sqaures: Variance [1]
out.savgh = zeros(1,size_glcm_3); % Sum average [1]
out.svarh = zeros(1,size_glcm_3); % Sum variance [1]
out.senth = zeros(1,size_glcm_3); % Sum entropy [1]
out.dvarh = zeros(1,size_glcm_3); % Difference variance [4]
%out.dvarh2 = zeros(1,size_glcm_3); % Difference variance [1]
out.denth = zeros(1,size_glcm_3); % Difference entropy [1]
out.inf1h = zeros(1,size_glcm_3); % Information measure of correlation1 [1]
out.inf2h = zeros(1,size_glcm_3); % Informaiton measure of correlation2 [1]
%out.mxcch = zeros(1,size_glcm_3);% maximal correlation coefficient [1]
%out.invdc = zeros(1,size_glcm_3);% Inverse difference (INV) is homom [3]
out.indnc = zeros(1,size_glcm_3); % Inverse difference normalized (INN) [3]
out.idmnc = zeros(1,size_glcm_3); % Inverse difference moment normalized [3]
% correlation with alternate definition of u and s
%out.corrm2 = zeros(1,size_glcm_3); % Correlation: matlab
%out.corrp2 = zeros(1,size_glcm_3); % Correlation: [1,2]
glcm_sum  = zeros(size_glcm_3,1);
glcm_mean = zeros(size_glcm_3,1);
glcm_var  = zeros(size_glcm_3,1);
% http://www.fp.ucalgary.ca/mhallbey/glcm_mean.htm confuses the range of 
% i and j used in calculating the means and standard deviations.
% As of now I am not sure if the range of i and j should be [1:Ng] or
% [0:Ng-1]. I am working on obtaining the values of mean and std that get
% the values of correlation that are provided by matlab.
u_x = zeros(size_glcm_3,1);
u_y = zeros(size_glcm_3,1);
s_x = zeros(size_glcm_3,1);
s_y = zeros(size_glcm_3,1);
% % alternate values of u and s
% u_x2 = zeros(size_glcm_3,1);
% u_y2 = zeros(size_glcm_3,1);
% s_x2 = zeros(size_glcm_3,1);
% s_y2 = zeros(size_glcm_3,1);
% checked p_x p_y p_xplusy p_xminusy
p_x = zeros(size_glcm_1,size_glcm_3); % Ng x #glcms[1]  
p_y = zeros(size_glcm_2,size_glcm_3); % Ng x #glcms[1]
p_xplusy = zeros((size_glcm_1*2 - 1),size_glcm_3); %[1]
p_xminusy = zeros((size_glcm_1),size_glcm_3); %[1]
% checked hxy hxy1 hxy2 hx hy
hxy  = zeros(size_glcm_3,1);
hxy1 = zeros(size_glcm_3,1);
hx   = zeros(size_glcm_3,1);
hy   = zeros(size_glcm_3,1);
hxy2 = zeros(size_glcm_3,1);
%Q    = zeros(size(glcm));
for k = 1:size_glcm_3 % number glcms
    glcm_sum(k) = sum(sum(glcm(:,:,k)));
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k); % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k)); % compute mean after norm
    glcm_var(k)  = (std2(glcm(:,:,k)))^2;
    
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            out.contr(k) = out.contr(k) + (abs(i - j))^2.*glcm(i,j,k);
            out.dissi(k) = out.dissi(k) + (abs(i - j)*glcm(i,j,k));
            out.energ(k) = out.energ(k) + (glcm(i,j,k).^2);
            out.entro(k) = out.entro(k) - (glcm(i,j,k)*log(glcm(i,j,k) + eps));
            out.homom(k) = out.homom(k) + (glcm(i,j,k)/( 1 + abs(i-j) ));
            out.homop(k) = out.homop(k) + (glcm(i,j,k)/( 1 + (i - j)^2));
            % [1] explains sum of squares variance with a mean value;
            % the exact definition for mean has not been provided in 
            % the reference: I use the mean of the entire normalized glcm 
            out.sosvh(k) = out.sosvh(k) + glcm(i,j,k)*((i - glcm_mean(k))^2);
            
            %out.invdc(k) = out.homom(k);
            out.indnc(k) = out.indnc(k) + (glcm(i,j,k)/( 1 + (abs(i-j)/size_glcm_1) ));
            out.idmnc(k) = out.idmnc(k) + (glcm(i,j,k)/( 1 + ((i - j)/size_glcm_1)^2));
            u_x(k)          = u_x(k) + (i)*glcm(i,j,k); % changed 10/26/08
            u_y(k)          = u_y(k) + (j)*glcm(i,j,k); % changed 10/26/08
            % code requires that Nx = Ny 
            % the values of the grey levels range from 1 to (Ng) 
        end
        
    end
    out.maxpr(k) = max(max(glcm(:,:,k)));
end
% glcms have been normalized:
% The contrast has been computed for each glcm in the 3D matrix
% (tested) gives similar results to the matlab function
for k = 1:size_glcm_3
    
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            p_x(i,k) = p_x(i,k) + glcm(i,j,k); 
            p_y(i,k) = p_y(i,k) + glcm(j,i,k); % taking i for j and j for i
            if (ismember((i + j),[2:2*size_glcm_1])) 
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcm(i,j,k);
            end
            if (ismember(abs(i-j),[0:(size_glcm_1-1)])) 
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) +...
                    glcm(i,j,k);
            end
        end
    end
    
%     % consider u_x and u_y and s_x and s_y as means and standard deviations
%     % of p_x and p_y
%     u_x2(k) = mean(p_x(:,k));
%     u_y2(k) = mean(p_y(:,k));
%     s_x2(k) = std(p_x(:,k));
%     s_y2(k) = std(p_y(:,k));
    
end
% marginal probabilities are now available [1]
% p_xminusy has +1 in index for matlab (no 0 index)
% computing sum average, sum variance and sum entropy:
for k = 1:(size_glcm_3)
    
    for i = 1:(2*(size_glcm_1)-1)
        out.savgh(k) = out.savgh(k) + (i+1)*p_xplusy(i,k);
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
        out.senth(k) = out.senth(k) - (p_xplusy(i,k)*log(p_xplusy(i,k) + eps));
    end
end
% compute sum variance with the help of sum entropy
for k = 1:(size_glcm_3)
    
    for i = 1:(2*(size_glcm_1)-1)
        out.svarh(k) = out.svarh(k) + (((i+1) - out.senth(k))^2)*p_xplusy(i,k);
        % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
    end
end
% compute difference variance, difference entropy, 
for k = 1:size_glcm_3
% out.dvarh2(k) = var(p_xminusy(:,k));
% but using the formula in 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
% we have for dvarh
    for i = 0:(size_glcm_1-1)
        out.denth(k) = out.denth(k) - (p_xminusy(i+1,k)*log(p_xminusy(i+1,k) + eps));
        out.dvarh(k) = out.dvarh(k) + (i^2)*p_xminusy(i+1,k);
    end
end
% compute information measure of correlation(1,2) [1]
for k = 1:size_glcm_3
    hxy(k) = out.entro(k);
    for i = 1:size_glcm_1
        
        for j = 1:size_glcm_2
            hxy1(k) = hxy1(k) - (glcm(i,j,k)*log(p_x(i,k)*p_y(j,k) + eps));
            hxy2(k) = hxy2(k) - (p_x(i,k)*p_y(j,k)*log(p_x(i,k)*p_y(j,k) + eps));
%             for Qind = 1:(size_glcm_1)
%                 Q(i,j,k) = Q(i,j,k) +...
%                     ( glcm(i,Qind,k)*glcm(j,Qind,k) / (p_x(i,k)*p_y(Qind,k)) ); 
%             end
        end
        hx(k) = hx(k) - (p_x(i,k)*log(p_x(i,k) + eps));
        hy(k) = hy(k) - (p_y(i,k)*log(p_y(i,k) + eps));
    end
    out.inf1h(k) = ( hxy(k) - hxy1(k) ) / ( max([hx(k),hy(k)]) );
    out.inf2h(k) = ( 1 - exp( -2*( hxy2(k) - hxy(k) ) ) )^0.5;
%     eig_Q(k,:)   = eig(Q(:,:,k));
%     sort_eig(k,:)= sort(eig_Q(k,:),'descend');
%     out.mxcch(k) = sort_eig(k,2)^0.5;
% The maximal correlation coefficient was not calculated due to
% computational instability 
% http://murphylab.web.cmu.edu/publications/boland/boland_node26.html
end
corm = zeros(size_glcm_3,1);
corp = zeros(size_glcm_3,1);
% using http://www.fp.ucalgary.ca/mhallbey/glcm_variance.htm for s_x s_y
for k = 1:size_glcm_3
    for i = 1:size_glcm_1
        for j = 1:size_glcm_2
            s_x(k)  = s_x(k)  + (((i) - u_x(k))^2)*glcm(i,j,k);
            s_y(k)  = s_y(k)  + (((j) - u_y(k))^2)*glcm(i,j,k);
            corp(k) = corp(k) + ((i)*(j)*glcm(i,j,k));
            corm(k) = corm(k) + (((i) - u_x(k))*((j) - u_y(k))*glcm(i,j,k));
            out.cprom(k) = out.cprom(k) + (((i + j - u_x(k) - u_y(k))^4)*...
                glcm(i,j,k));
            out.cshad(k) = out.cshad(k) + (((i + j - u_x(k) - u_y(k))^3)*...
                glcm(i,j,k));
        end
    end
    % using http://www.fp.ucalgary.ca/mhallbey/glcm_variance.htm for s_x
    % s_y : This solves the difference in value of correlation and might be
    % the right value of standard deviations required 
    % According to this website there is a typo in [2] which provides
    % values of variance instead of the standard deviation hence a square
    % root is required as done below:
    s_x(k) = s_x(k) ^ 0.5;
    s_y(k) = s_y(k) ^ 0.5;
    out.autoc(k) = corp(k);
    out.corrp(k) = (corp(k) - u_x(k)*u_y(k))/(s_x(k)*s_y(k));
    out.corrm(k) = corm(k) / (s_x(k)*s_y(k));
%     % alternate values of u and s
%     out.corrp2(k) = (corp(k) - u_x2(k)*u_y2(k))/(s_x2(k)*s_y2(k));
%     out.corrm2(k) = corm(k) / (s_x2(k)*s_y2(k));
end
% Here the formula in the paper out.corrp and the formula in matlab
% out.corrm are equivalent as confirmed by the similar results obtained
% % The papers have a slightly different formular for Contrast
% % I have tested here to find this formula in the papers provides the 
% % same results as the formula provided by the matlab function for 
% % Contrast (Hence this part has been commented)
% out.contrp = zeros(size_glcm_3,1);
% contp = 0;
% Ng = size_glcm_1;
% for k = 1:size_glcm_3
%     for n = 0:(Ng-1)
%         for i = 1:Ng
%             for j = 1:Ng
%                 if (abs(i-j) == n)
%                     contp = contp + glcm(i,j,k);
%                 end
%             end
%         end
%         out.contrp(k) = out.contrp(k) + n^2*contp;
%         contp = 0;
%     end
%     
% end
%       GLCM Features (Soh, 1999; Haralick, 1973; Clausi 2002)
%           f1. Uniformity / Energy / Angular Second Moment (done)
%           f2. Entropy (done)
%           f3. Dissimilarity (done)
%           f4. Contrast / Inertia (done)
%           f5. Inverse difference    
%           f6. correlation
%           f7. Homogeneity / Inverse difference moment
%           f8. Autocorrelation
%           f9. Cluster Shade
%          f10. Cluster Prominence
%          f11. Maximum probability
%          f12. Sum of Squares
%          f13. Sum Average
%          f14. Sum Variance
%          f15. Sum Entropy
%          f16. Difference variance
%          f17. Difference entropy
%          f18. Information measures of correlation (1)
%          f19. Information measures of correlation (2)
%          f20. Maximal correlation coefficient
%          f21. Inverse difference normalized (INN)
%          f22. Inverse difference moment normalized (IDN)
end


function [x] = haralickTextureFeatures(coOcMat, xFeatures)
%https://ch.mathworks.com/matlabcentral/fileexchange/58769-haralicktexturefeatures
%Calculates all Haralick Features.
%
% Function call:
%   [x] = haralickTextureFeatures(coOcMat) calculates all 14 Haralick
%   Features
%	[x] = haralickTextureFeatures(coOcMat, xFeatures) calculates the
%	Haralick Features specified by xFeatures, the rest will be return as 0.
%	Use this for better legacy if you do not need all Haralick Features.
%
%
% Source:           http://haralick.org/journals/TexturalFeatures.pdf
%
%
% input:
%   'coOcMat'       Co-Occurence-Matrix,  which must be a [nxm] matrix,
%                   see matlab documentation glcm
%   'xFeatures'     (optional) - Feature(s), which should be calculated
%
% output:           
%   'x' - [vector with the following feature(s):
%               x(1)  Angular Second Moment (Energy) [checked]
%               x(2)  Contrast [checked]
%               x(3)  Correlation [checked]
%               x(4)  Variance [checked]
%               x(5)  Inverse Difference Moment (Homogeneity) [checked]
%               x(6)  Sum Average [checked]
%               x(7)  Sum Variance [approxemitly (cut out zeros)]
%               x(8)  Sum Entropy [checked]
%               x(9)  Entropy [cut out zeros]
%               x(10) Difference Variance [approxemitly]
%               x(11) Difference Entropy [checked]
%               x(12) Information Measure of Correlation I [checked]
%               x(13) Information Measure of Correlation II [approxemitly]
%               x(14) Maximal Correlation Coefficient [validated, no reference]
%
%
%
%   Example
%   ---------      
%   %Load an image of Matlab
%   I = imread('circuit.tif');
% 
%	%Get the co occurence matrix (in Matlab called GLCM: Gray Level Co 
%   %Occurence Matrix)
%   glcm = graycomatrix(I, 'offset', [0 1], 'Symmetric', true);
% 
%	%calculate feature 4 (Variance), 5 (Inverse Difference Moment) and 6 
%   %(Sum Average)
%   xFeatures = 4:6;
%   x = haralickTextureFeatures(glcm, 4:6);
%
%   %Get only the features you want
%   x = x( xFeatures )    
%
%
% Notes:        If x14 Maximal Correlation Coefficient is complex then the
%               magnitude of MCC will be calculate.
%               See the haralick paper to understand the code.
%
% Info:         ver 1.1
%               - coOcMat will be checked if it is 2-dimensional
%               - Example code added
%               - more documentation
%               - fixed if-end polling (thanks to Ihsan Yassin)
%
%               ver 1.2
%               - fixed and validated Haralick Feature 14, still no
%               reference
%               - fixed problem with empty eigenvec matrix (thanks to Hafiz
%               Muhammad Arslan)
%
%               ver 1.3
%               - fixed feature 8 Sum Variance (thanks to Lingxuan Kong)
%
% Author:       Rune Monzel, runemonzel(at)gmail.com
%
% See also graycomatrix, graycoprops.
% check input
if nargin == 1
    xFeatures = 1 : 14;
end
% check coOcMat for dimensions:
if ~(ismatrix(coOcMat))
    error(['\coOcMatInput must be a two dimensional matrix, '...
        'dimensional was %s.',ndims(coOcMat)']);
end
% initialize x
x = zeros(14,1);
% normalize glcm
coOcMat = coOcMat./sum(coOcMat(:));
%% Some pre-calculation:
% columns and rows
if sum(xFeatures == 2) == 1 | ... % Contrast
        sum(xFeatures == 3) == 1 | ... % Correlation
        sum(xFeatures == 4) == 1 | ... % Variance
        sum(xFeatures == 5) == 1 | ... % Inverse Difference Moment
        sum(xFeatures == 6) == 1 | ... % Sum Average
        sum(xFeatures == 7) == 1 | ... % Sum Variance
        sum(xFeatures == 8) == 1 | ... % Sum Entropy
        sum(xFeatures == 10) == 1 | ...% Difference Variance
        sum(xFeatures == 11) == 1 | ...% Difference Entropy
        sum(xFeatures == 14) == 1 % Maximal Correlation Coefficient
    sizecoOcMat = size(coOcMat);
    [col,row] = meshgrid(1:sizecoOcMat(1),1:sizecoOcMat(2));
end
% average and standarddeviation
if sum(xFeatures == 3) == 1 | ... % correlation
        sum(xFeatures == 10) == 1 % difference variance
    
    
    rowMean =  sum( row(:).*coOcMat(:) );
    colMean = sum( col(:).*coOcMat(:) );
    rowStd = sqrt( sum( (row(:)-rowMean).^2 .* coOcMat(:) ) );
    colStd = sqrt( sum( (col(:)-colMean).^2 .* coOcMat(:) ) );
end
% sum of rows p_y(i) und sum of columns p_x(j)
if sum(xFeatures == 12) == 1 |...% Information Measures of Correlation I
        sum(xFeatures == 13) == 1|... % Information Measures of Correlation II
        sum(xFeatures == 14) == 1 % Maximal Correlation Coefficient
    
    rowCoOcMat = sum(coOcMat,2); %sum of rows p_y(i)
    colCoOcMat = sum(coOcMat); %sum of columns p_x(i)
end
% p_x+y
if sum(xFeatures == 6)==1 |... % Sum Average
        sum(xFeatures == 7)==1 |... % Sum Variance
        sum(xFeatures == 8)==1 % Sum Entropy
    
    start = -(sizecoOcMat(1) -1);
    stop = sizecoOcMat(1) -1;
    
    % Rotate Matrix 90??
    coOcMat90 = rot90(coOcMat);
    
    % Initilisiere p_x+y
    p_XplusY = zeros((2*sizecoOcMat(1))-1,1);
    
    k = 1;
    for index = start : stop
        p_XplusY(k) = sum( diag(coOcMat90,index) );
        k = k + 1;
    end
end
% Initialize  p_x-y
if sum(xFeatures == 10)==1 |... % Difference Variance
        sum(xFeatures == 11)==1 % Difference Entropy
    
    start = 1;
    stop = sizecoOcMat(1)-1;
    
    % Initialize p_XminusY
    p_XminusY = zeros(sizecoOcMat(1),1);
    p_XminusY(1) = sum (diag(coOcMat,0) );
    
    k = 2;
    for index = start : stop
        p_XminusY(k) = sum( [diag(coOcMat,index);
            diag(coOcMat,-index)] );
        k = k + 1;
    end
end
%% Haralick Feature Calculations
for f = xFeatures
    switch f
        case 1 % Energy (Angular Second Moment)
            x(1) = sum( coOcMat(:).^2 );
            
        case 2  % Contrast
            matrix = ( abs(row - col).^2 ) .* coOcMat;
            x(2) = sum( matrix(:) );
            
        case 3  % Correlation
            zaehler = sum ((row(:) - rowMean) .*...
                (col(:) - colMean) .*  coOcMat(:));
            denominator = rowStd * colStd;
            x(3) = zaehler/denominator;
            
        case 4 % Variance
            x(4) = sum( (row(:)-mean(coOcMat(:))).^2 .*coOcMat(:) );
            
        case 5 % Inverse Difference Moment
            x(5) = sum( coOcMat(:) ./ ( 1+ (row(:)-col(:)).^2 ) );
            
        case 6 % Sum Average
            x(6) = sum( (2:(2*sizecoOcMat(1)))' .* p_XplusY );
            
        case 7 % Sum Variance
            x(8) = - sum( p_XplusY(p_XplusY~=0) .* ...
                log(p_XplusY(p_XplusY~=0)) );
            
            x(7) = sum( ((2:(2*sizecoOcMat(1)))' -...
                x(8)).^2 .* p_XplusY  );
            
        case 8 % Sum Entropy
            if ~x(8) % if it is not calculate in case 7
                x(8) = - sum( p_XplusY(p_XplusY~=0) .*...
                    log(p_XplusY(p_XplusY~=0)) );
            end
            
        case 9 % Entropy
            x(9) = - sum( coOcMat(coOcMat~=0) .*...
                log2(coOcMat(coOcMat~=0)) );
            
        case 10 % Difference Variance
            x(10) = sum( ((0:sizecoOcMat(1)-1)' -...
                mean(p_XminusY)).^2 .* p_XminusY);
            
        case 11 % Difference Entropy
            x(11) = - sum( p_XminusY(p_XminusY~=0) .*...
                log(p_XminusY(p_XminusY~=0)) );
            
        case 12 % Information Measures of Correlation I
            
            x(9) = - sum( coOcMat(coOcMat~=0) .*...
                log2(coOcMat(coOcMat~=0)) );
            
            % Cuto out all zeros:
            logrc  = log2( rowCoOcMat*colCoOcMat ); % 256x1 * 1x256
            %Matrixmultiplication
            logrc(logrc == -Inf) = 0; % cut out Inf
            HXY1 = - sum( coOcMat(:).* logrc(:) ); %product of elements
            % between co-occurence-matrix and the logarithmetic matrix
            numerator = x(9) - HXY1;
            
            % calculate off HX, Entropy of sum of columns
            logc = log2(colCoOcMat);
            logc(logc==-Inf) = 0;
            HX = - sum( colCoOcMat .* logc );
            
            % calculate off HY, Entropy of sum of columns
            logr = log2( rowCoOcMat );
            logr(logr==-Inf) = 0;
            HY = - sum( rowCoOcMat .* logr );
            
            % max value
            denominator = max([HX HY]);
            x(12) = numerator / denominator;
            
        case 13 % Information Measures of Correlation II
            if x(9)
                x(9) = - sum( coOcMat(coOcMat~=0) .*...
                    log2(coOcMat(coOcMat~=0)) );
            end
            logrc  = log2( rowCoOcMat*colCoOcMat ); % 256x1 * 1x256
            %Matrixmultiplication
            logrc(logrc == -Inf) = 0;
            HXY2 = - sum( sum( (rowCoOcMat * colCoOcMat) .* logrc ));
            x(13) =  (  ( 1 - exp(-2*(HXY2 - x(9))) )  ).^(1/2);
            
        case 14 % Maximal Correlation Coefficient
            
            % Initialise Q
            Q = zeros(sizecoOcMat(1),sizecoOcMat(2));
            
            % % Slow version
            %for i = 1 : length(coOcMat(:,1))
            %    for j = 1 : length(coOcMat(1,:))
            %        Q(i,j) = sum(( coOcMat(i,:) .* coOcMat(j,:) )...
            %            ./ ( rowCoOcMat(i) .* colCoOcMat ),'omitnan');
            %    end
            %end
            
            % Fast version
            for i = 1 : sizecoOcMat(2)
                    Q(i,:) = sum( ...
                        (repmat(coOcMat(i,:),sizecoOcMat(1),1) .* coOcMat ) ./ ...
                        repmat( rowCoOcMat(i) .* colCoOcMat, sizecoOcMat(1),1),...
                        2,'omitnan');
            end
            
            % cut out nans
            Q(isnan(Q)) = 0;
            
            eigenvec = eig(Q);
            
            % Find largest eigenvec and delete
            eigenvec(eigenvec==max(eigenvec))=[];
            
            % If eigenvec is a zero matrix, no maximum can be found, that
            % would mean, that the second largest eigenvec of course is
            % also 0
            if ~any(eigenvec) || isempty(eigenvec)
                x(14) = sqrt(0);
                continue;
            end
            
            % Sqrt of second largest eigenvec
            x(14) = sqrt( max(eigenvec) );
            
            % calculate magnitude of Maximal Correlation Coefficient
            if imag(x(14))
                x(14) = abs(x(14));
            end
            
    end
end
end