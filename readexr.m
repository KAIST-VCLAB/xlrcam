function [RGB, fileinfo] = readexr(filename)
    % this function reads only rgb
    mapObj = exrreadchannels(filename);
 
   % sRGB
    RGB(:,:,1) = mapObj('R');
    RGB(:,:,2) = mapObj('G');
    RGB(:,:,3) = mapObj('B');
end