% converts grayscale image to BW mask
% for use with maskingGUI.m, maskingGUI.fig
% Mask methods are optimized for fluorescence images unless otherwise stated
% Written by Sophie Skanchy in version R2019b

% Masks list:
% Mask 1: Amir's Method
% Mask 2: Sobel Edge
% Mask 3: Canny + log Edge
% Mask 4: log Edge
% Mask 5: Median Filtering and imbinarize
% Mask 6: Majority of masks 1-5 agree
% Mask 7: at least 2 of masks 1-5 agree


function mask = im2mask(I, method)
% img is original grayscale image
% method is chosen masking method

mask_size = size(I);
mask = cell(1, 10);

%% Mask 1: Amir's Method
level=ceil(graythresh(I)*100)/100;
BW=im2bw(I,level);

% cleaning the BW image
%BW = bwmorph(BW, 'clean');
BW2 = imfill(BW, 'holes');
BW3 = bwmorph(BW2,'majority', 2);
BW4 = bwareaopen(BW3, 5);
BW5 = imclearborder(BW4);
BW6 = bwmorph(BW5,'thicken',1);

mask{1} = BW6;


%% Mask 2: Edge Detection - Sobel
edgeBW = edge(I, 'Sobel');
%close gaps
radius = 2; num = 4;
se = strel('disk', radius, num);
edgeBWc = imclose(edgeBW,se);

%fill interior
edgeBWfill = imfill(edgeBWc,'holes');

BW2 = imfill(edgeBWfill, 'holes');
BW3 = bwmorph(BW2,'majority', 2);
BW4 = bwareaopen(BW3, 5);
BW5 = imclearborder(BW4);
BW6 = bwmorph(BW5,'thicken',1);

mask{2} = BW6;


%% Mask 3: Edge Detection - Canny + log, no strel
% edge detection
edgeCanny = edge(I, 'Canny');
edgeLog = edge(I, 'log');

% fill and remove noise
fillCanny = imfill(edgeCanny, 'holes');
majCanny = bwmorph(fillCanny, 'majority');

fillLog = imfill(edgeLog, 'holes');
majLog = bwmorph(fillLog, 'majority');

% sum majCanny and majLog, remove noise
majCannyLog = (majCanny + majLog) > 0;
majCannyLog = bwmorph(majCannyLog,'majority', 3);

mask{3} = majCannyLog;


%% Mask 4: Edge Detection - log, no strel, minimal cleaning
edgeLog = edge(I, 'log');
fillLog = imfill(edgeLog, 'holes');
majLog = bwmorph(fillLog, 'majority', 2);

mask{4} = majLog;


%% Mask 5: Median Filtering and Graythresh
binarized = imbinarize(I, 'adaptive');
filt = medfilt2(binarized);
cleaned = bwmorph(filt,'majority', 2);

mask{5} = cleaned;


%% Mask 6: Majority of masks 1-5
sum_masks = mask{1} + mask{2} + mask{3} + mask{4} + mask{5};
average = sum_masks >= 3; % at least 3 of 5 masks agree
clean = bwmorph(average,'majority',2);

mask{6} = clean;


%% Mask 7: At least 2 of masks 1-5
sum_masks = mask{1} + mask{2} + mask{3} + mask{4} + mask{5};
average = sum_masks >= 2; % at least 2 of 5 masks agree
clean = bwmorph(average,'majority',2);

mask{7} = clean;




%% select mask
mask = mask{method};