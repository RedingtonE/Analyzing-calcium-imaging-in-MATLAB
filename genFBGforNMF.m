function [F, BG, mBG] = genFBGforNMF(masks, vidname)
%The goal of this function is to isolate the fluorescence activity of a
%set of neurons in a MxNxP matrix of regions of interest from a h5 video of calcium activity.
%Input: 
%masks = a MxNxP matrix of regions of interest from P neurons
%vidname = the name of the h5 video to calculate the fluorescence activity.
%The video itself is an MxNxT matrix with T frames
%Output:
%F = a PxT matrix of the average fluorescence of the ROI
%BG = a PxT matrix of the average fluorescence of the ROI
%mBG = a PxT matrix of the median background fluorescence of the ROI
%from
%Load in the video
   video = h5read(vidname, '/regvideo');
%Create a single image that is a composite of all the neurons. This is then
%used to identify all non-neuron pixels in the field of view that will be
%used as candidates for calculating the background of the video. 
   allMask = sum(masks, 3);
   allMask(allMask > 0) = 1;
   allMask = ~allMask;
   mallMask = ones(size(allMask));
   [d1 d2 d3] = size(video);
%Make a (MxN)xT version of the video that will be used later in the code. I
%do this to speed up calculating the average fluorescence on a frame by
%frame basis. Rather than running a for loop for every frame for every
%neuron, I can instead take advantage of how MATLAB is optimized for matrix
%functions and calculate the average through all of T. This speeds up the
%calculation significantly. 
   linvid = reshape(video, [d1*d2 d3]);
   [xx yy] = meshgrid(1:d1,1:d2);
   BG = zeros(d3, size(masks,3));
   mBG = zeros(d3, size(masks, 3));
   F = zeros(d3, size(masks, 3));
   for fn = 1:size(masks, 3)
       %Identify the size of the neuron. This will be used to determine the
       %number of pixels used in the background
        stats = regionprops(logical(masks(:, :, fn)),'MajorAxisLength','MinorAxisLength');
        nsize = (stats.MajorAxisLength+stats.MinorAxisLength)/2;
        radi = nsize/2;
        radi2 = nsize/2;
        
        %Identify the COM of the region of interest and initialize the
        %BGmask that will be used to calculate the background fluorescence
        %surrounding that neuron.
        [X Y] = findCOM(masks(:, :, fn));
        BGmask = zeros(size(masks(:, :, fn)));
        BGsize = length(find(BGmask));
        mBGmask = zeros(size(masks(:, :, fn)));
        %Grow the size of the background mask until it has at least the
        %desired number of pixels
        while BGsize <= length(find(masks(:, :, fn)))*1
        x1 = max([1 round(X-radi)]);
        x2 = min([d1 round(X+radi)]);
        y1 = max([1 round(Y - radi)]);
        y2 = min([d2 round(Y+radi)]);
        
        circleout = double((yy-X).^2 + (xx-Y).^2 < radi^2)';
        BGmask = circleout.*allMask;
        BGsize = length(find(BGmask));
        radi = radi + 0.5;
        end
        
        BGind = find(BGmask);
        maskind = find(masks(:, :, fn));

        %Take advantage of the linear video to easily use the indexes of
        %the background mask and neural mask (region of interest)
        BGtotal = sum(linvid(BGind, :),1);
        mBG(:, fn) = median(linvid(BGind, :), 1);
        Ftot = sum(linvid(maskind, :),1);
        BG(:, fn) = BGtotal/BGsize;
        F(:, fn) = Ftot/length(find(masks(:, :, fn)));
   end
end
        
        
