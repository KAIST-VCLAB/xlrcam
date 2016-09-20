%=========================================================================%
% ``Modeling Human Color Perception under Extended Luminance Levels''
% ACM Transactions on Graphics, presented in SIGGRAPH 2009.
% 
% This is a high-dynamic-range color appearance model implemented in Matlab
% ver. 1.5 (released in 22/06/2009; last update in 01/19/2016)
%
% [Reference]
% @Article{KimWeyKautz:2009:SIG,
%   author  = {Min H. Kim and Tim Weyrich and Jan Kautz},
%   title   = {Modeling Human Color Perception
%              under Extended Luminance Levels},
%   journal = {ACM Transactions on Graphics (Proc. SIGGRAPH 2009)},
%   year    = {2009},
%   volume  = {28},
%   number  = {3},
%   pages   = {27:1--9},
%   doi     = "1531326.1531333",
%   URL     = "http://dl.acm.org/citation.cfm?doid=1531326.1531333"
% }
%=========================================================================%
% [How to use this code for a single HDR image]
% Just call "xlrcam.m" for general purposes
% xlrcam                
%
% [Usages with characterized HDR images in absolute scales]
% 1. relative chromaticity mapping space
% xlrcam('JCh','noCAT') % JCh: for relative HDR images (general)
% xlrcam('JCh','CAT')   %                "             plus CAT
% xlrcam('JCh','CAT',hdrscaling) % applying a scaling constant to HDR images
%
% 2. absolute chromaticity mapping space
% xlrcam('JMh','noCAT') % JMh: for absolute HDR images (calibrated)
% xlrcam('JMh','CAT')   %                "             plus CAT
% xlrcam('JMh','CAT',hdrscaling) % applying a scaling constant to HDR images
%
% Note that CAT is implemented with automatic white-point estimation,
% which take the brightest pixel value in the image. Thus its estimate 
% of the image reference white could be inaccurate.
%=========================================================================%
% [Example code for HDR video frames (trick to avoid flickering)]
%     for i = 0:( num_images - 1 )
%         frame_idx = start_frame_index + i;
%         % Read images
%         filename_to_read = sprintf( '%s%d%s', prefix_previous_job, frame_idx, file_extension_HDR );
%         fullfile_to_read = fullfile( out_directory, dir_previous_job, filename_to_read );
%         [ HDR_image, fileinfo] = readhdr( fullfile_to_read );
%         fprintf( '[HDR-TM] Processing %d th frame...( %d/%d )\n', frame_idx, i + 1, num_images );
%         if i == 0 % we estimate adaptive parameters 
%             [ HDR_TM_image, La, inXYZwhite, XYZ2ma ] = hlcamrep( HDR_image, outXYZw, outLa, outmedia, CSpace, CAT );
%         else
%             [ HDR_TM_image, La, inXYZwhite, XYZ2ma ] = hlcamrep( HDR_image, outXYZw, outLa, outmedia, CSpace, CAT, La, inXYZwhite, XYZ2ma );
%         end
%         if save_TM_HDR
%             img_to_save = HDR_TM_image;
%             file_name_to_save = sprintf( '%s%d%s', prefix_TM_HDR, frame_idx, file_extension );
%             fullfile_to_save = fullfile( out_directory, dir_TM_HDR, file_name_to_save );
%             imwrite( img_to_save, fullfile_to_save );
%         end
%     end
%=========================================================================%
% Version history
% 1.0: 22/06/2009 : first release for SIGGRAPH 2009
% 1.1: 17/01/2015 : update in order to make this code handle rendered images (with zeros) and HDR video tone-mapping (to avoid flickering)
% 1.2: 16/10/2015 : suport for OpenEXR
% 1.3: 20/10/2015 : library link fixed for OpenEXR
% 1.4: 01/18/2016 : adding scale function
% 1.5: 01/18/2016 : adding media option
%=========================================================================%
%     Copyright (c) 2009-16, Min H. Kim
%     All rights reserved.
% 
%     Min H. Kim has developed "xlrcam.m" and related documentation 
%     (the "Software"); confidential use in source form of the Software, 
%     without modification, is permitted provided that the following 
%     conditions are met:
%     1. Neither the name of the copyright holder nor the names of any 
%     contributors may be used to endorse or promote products derived from 
%     the Software without specific prior written permission.
%     2. The use of the software is for Non-Commercial Purposes only. As 
%     used in this Agreement, "Non-Commercial Purpose" means for the 
%     purpose of education or research in a non-commercial organisation 
%     only. "Non-Commercial Purpose" excludes, without limitation, any use 
%     of the Software for, as part of, or in any way in connection with a 
%     product (including software) or service which is sold, offered for 
%     sale, licensed, leased, published, loaned or rented. If you require 
%     a license for a use excluded by this agreement, 
%     please email [minhkim@kaist.ac.kr].
%
%     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
%     ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
%     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
%     A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR 
%     CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
%     EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
%     PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
%     PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY 
%     OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
%     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
%     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%=========================================================================%
% OpenEXR 1.7.0 compiled as MEX with courtesy of:
% Edgar Velazquez-Armendariz (eva5@cs.cornell.edu)
% Based on the originals by Jinwei Gu (jwgu@cs.columbia.edu)
%=========================================================================%
    
function xlrcam(CSpace,CAT,hdrscaling,outmedia)

    %=====================================================================%
    % Default Parameters for an output target medium:
    outXYZw = D(65)*2.5;        % sRGB, 250cd/sqm (Apple Cinema HD display)
    outLa = 0.10 * outXYZw(2);  % 25 cd/sqm (10% of peak lumin.) ~= dim
    % outmedia options: hdr, trans, lcd, crt, paper
    % outmedia defaults: high quality LCD (between CRT and HDR)
    %=====================================================================%    
    if nargin<1
       CSpace = 'JCh';
       CAT = 'noCAT';
       hdrscaling = 1.0;
       outmedia = 'trans';
    elseif nargin<2
       CAT = 'noCAT';
       hdrscaling = 1.0;
       outmedia = 'trans';
    elseif nargin<3
       hdrscaling = 1.0;
       outmedia = 'trans';
    elseif nargin<4
       outmedia = 'trans';
    end
    % file dialogue
    [file_name,file_path] = uigetfile({'*.hdr;*.exr','HDR'},'Choose HDR images','MultiSelect', 'on');
    if iscell(file_name)==0;
        if file_name==0; 
            disp(sprintf('No file selected'));
            return; 
        end;
    end
    source_images = strcat(file_path,file_name);
    filenames = char(source_images);
    filenames = sortrows(filenames);

    global filename;
    whole = size(filenames,1);
    for j=1 : whole
        filename = strtrim(filenames(j,:));
        disp(sprintf('==================================='));
        disp(sprintf('Filename   : %s', filename ));
        if ~isempty(findstr(filename, '.hdr'))
            inimage = hdrscaling.*readhdr(filename);
        elseif ~isempty(findstr(filename, '.exr'))
            inimage = hdrscaling.*readexr(filename);
        else
            return;
        end
%        inimage = 50.*inimage;%scaling luminance
        outimage = hlcamrep(inimage,outXYZw,outLa,outmedia,CSpace,CAT);
%        outimage = stretch_hist(outimage); % histgram streching
        fig = figure; imshow(outimage);
        if strcmp(CSpace,'JCh')
            if strcmp(CAT,'CAT')
                set(fig,'Name', 'JCh mode with CAT');
                imwrite(outimage, strcat(filename,'_xlrcam(JCh,CAT).png'));
            elseif  strcmp(CAT,'noCAT')
                set(fig,'Name', 'JCh mode with no CAT');
                imwrite(outimage, strcat(filename,'_xlrcam(JCh,noCAT).png'));
            end
        elseif strcmp(CSpace,'JMh')
            if strcmp(CAT,'CAT')
                set(fig,'Name', 'JMh mode with CAT');
                imwrite(outimage, strcat(filename,'_xlrcam(JMh,CAT).png'));
            elseif  strcmp(CAT,'noCAT')
                set(fig,'Name', 'JMh mode with no CAT');
                imwrite(outimage, strcat(filename,'_xlrcam(JMh,noCAT).png'));
            end
        end
        clear('inimage');
        clear('outimage');
        % close all; % uncomment for video processing
    end
end

function [ outimage, La, inXYZwhite, XYZ2max ] = hlcamrep(inimage,outXYZw,outLa,outmedia,CSpace,CAT, La, inXYZwhite, XYZ2max)
    % color transform
    MPCS = [0.4361,0.3851,0.1431;
        0.2225,0.7169,0.0606;
        0.0139,0.0971,0.7141;];

    sizeimg = size(inimage);
    inXYZ = (MPCS*(m2v(inimage))')';
    ndim = size(inXYZ,1);
    
    if nargin < 7
        inXYZtmp = m2v(imresize(v2m(inXYZ,size(inimage)),0.3,'bicubic'));
        [CC, NN] = max(inXYZtmp(:,2));
        inXYZwhite = inXYZtmp(NN,:);

        FLT_MAX = 1.0E+30;
        FLT_MIN = 1./1.0E+30;
        Y_max = -FLT_MAX;
        Y_min = FLT_MAX;
        for n = 1:size(inXYZtmp,1)
            temp = inXYZtmp(n,2);
            if (temp>Y_max) 
                Y_max = temp;
            end
            if ((temp > 0)&&(temp<Y_min)) 
                Y_min = temp;
            end
        end
    else
        Y_max = inXYZwhite(2);
    end
    

    disp(sprintf('==================================='));
    if strcmp(CAT,'CAT')
        disp(sprintf('With Chromatic Adaptation'));
    elseif strcmp(CAT,'noCAT')
        disp(sprintf('Without Chromatic Adaptation'));
    else
        disp(sprintf('Input Error!'));
    end
    
    if nargin < 7
        % La bug fix (2015-01-17): exclude zero levels in calculating La (for rendered images)
        Lworld = inXYZ(:,2);
        Lworld_nozeros = Lworld(Lworld~=0);
        Lworld_nozeros = Lworld_nozeros + FLT_MIN;
        log_Lw_nozeros = log(Lworld_nozeros);
        La = exp(mean(log_Lw_nozeros)); 
    end
    
%     Lworld = inXYZ(:,2)+FLT_MIN;
%     log_Lw = log(Lworld);
%     La = exp(mean(log_Lw));

    XYZ(1,:) = inXYZwhite;
    XYZ(3:ndim+2,:) = inXYZ;
    XYZ(2,2) = La;
    [J,C,h, Q,M,H, ac,bc,s]=xlrcamf(XYZ, 'hdr', CAT); % forward
    [XYZ2] = xlrcami(J, C, h, Q, M, outXYZw, outLa, outmedia, CAT,sizeimg, CSpace); % inverse

    if nargin < 7
        XYZ2tmp = m2v(imresize(v2m(XYZ2,sizeimg),0.3,'bicubic')); % ignoring noise, estimating white
        [CC, NN] = max(XYZ2tmp(:,2));
        XYZ2max = XYZ2tmp(NN,:); 
    end

    XYZ2n = 100.*XYZ2./XYZ2max(2);
    srgb3 = xyzd502srgb(XYZ2n);
    outimage  = uint8(v2m(srgb3,sizeimg));
end

function RGB=xyzd502srgb(XYZ0)
    M=[  3.1336,-1.6168,-0.4907;
        -0.9787, 1.9161, 0.0335;
         0.0721,-0.2291, 1.4054;];

    if ischar(XYZ0)
       XYZ=dlmread(XYZ0,'\t');
    elseif isnumeric(XYZ0)
       XYZ=XYZ0;
    else
       error('No valid input data')
    end

    if ndims(XYZ0)==3
        row=size(XYZ0,1);
        col=size(XYZ0,2);
        XYZ = reshape(XYZ0, row*col,3);
    end

    sRGB=(M*(XYZ./100)')';

    sR=sRGB(:,1);sG=sRGB(:,2);sB=sRGB(:,3);
    sR(sR>1)=1;sG(sG>1)=1;sB(sB>1)=1;
    sR(sR<0)=0;sG(sG<0)=0;sB(sB<0)=0;

    j=find(sR<=0.00304);
    k=find(sG<=0.00304);
    l=find(sB<=0.00304);

    g=1/2.4;
    R=(1.055*sR.^g-0.055)*255;
    G=(1.055*sG.^g-0.055)*255;
    B=(1.055*sB.^g-0.055)*255;    
    R(j)=(sR(j)*12.92)*255;
    G(k)=(sG(k)*12.92)*255;
    B(l)=(sB(l)*12.92)*255;

    R(R>255)=255;G(G>255)=255;B(B>255)=255;
    R(R<0)=0;G(G<0)=0;B(B<0)=0;
    RGB=[R,G,B];
    RGB = uint8(RGB);

    if ndims(XYZ0)==3
        RGB=reshape(RGB, row,col,3);
    end
end

function M=v2m(V,array_size)
    r=array_size(1);
    c=array_size(2);
    [vr,vc]=size(V);

    if vr~=r*c
       disp('Error in array size')
       return
    end
    M=zeros(r,c,vc);
    
    for i=1:vc
       Vm=V(:,i);
       Mt=reshape(Vm,c,r);
       M(:,:,i)=transpose(Mt);
    end
end

function v = m2v(data)
    LMN=data;
    [r,c,n]=size(LMN);
    v=zeros(r*c,n);

    for i=1:n
       L=LMN(:,:,i);
       Lt=L';
       v(:,i)=Lt(:);
    end
end

function XYZ=D(K)
       A=[109.850,100,35.585];
       C=[98.074,100,118.232];
       D50=[96.422,100,82.521];
       D55=[95.682,100,92.149];
       D65=[95.047,100,108.883];
       D75=[94.96,100,122.62];
       allXYZ=[A;C;D50;D55;D65;D75];
       xyz=num2str(allXYZ,'%8.2f');
       labels=[' A    ';' C    ';' D50  ';' D55  ';' D65  ';' D75  '];
       all=[labels,xyz];
    if nargin>0
       switch K
       case 'cie_a';XYZ=A;
       case 'cie_c';XYZ=C;
       case 50;XYZ=D50;
       case 55;XYZ=D55;
       case 60;XYZ=D60;
       case 65;XYZ=D65;
       case 75;XYZ=D75;
       end
    else
       disp(all);
       if nargout>0;XYZ=allXYZ;end
    end   
end

% Forward Model
function [J,C,h, Q,M,H, ac,bc,s]=xlrcamf(xyz, media, CAT);
% transform matrices
MCAT02=[0.7328,0.4296,-0.1624;
-0.7036,1.6975,0.0061;
0.0030,0.0136,0.9834];
MH=[0.38971,0.68898,-0.07868;
-.22981,1.1834,0.04641;
0,0,1];
MCAT02i=[1.096124,-0.278869,0.182745;
0.454369,0.473533,0.072098;
-0.009628,-0.005698,1.015326];
if nargin<2
   media = 'hdr';
   CAT = 'CAT';
elseif nargin<3
   CAT = 'CAT';
end

if strcmp(media,'hdr')
    mda = 1.0000;
elseif strcmp(media,'trans')
    mda = 1.2175;
elseif strcmp(media,'lcd')
    mda = 1.3374;
elseif strcmp(media,'crt')
    mda = 1.4572;
elseif strcmp(media,'paper')
    mda = 1.7526;
end

ndim=size(xyz);
XYZw = xyz(1,:);  Lw = xyz(1,2);
La = xyz(2,2);
XYZ = xyz(3:ndim(1),:);

% if you want to calculate the human color perception of the experimental
% data (a disk stimulus surrounded by background), please use this La calculation. 
% La = 0.88 * La + 0.08* mean(XYZ(:,2)) + 0.04 * Lw; 
%
% Otherwise for general tonemapping applications with HDR images
La = La;

MHM=MH*MCAT02i;

RGBw =(MCAT02*XYZw')';
RGB =(MCAT02*XYZ')';

if strcmp(CAT,'CAT')
    RGBwp(1) = RGBw(1)/(RGBw(1)/Lw);
    RGBwp(2) = RGBw(2)/(RGBw(2)/Lw);
    RGBwp(3) = RGBw(3)/(RGBw(3)/Lw);
    RGBp(:,1) = RGB(:,1)./(RGBw(1)/Lw);
    RGBp(:,2) = RGB(:,2)./(RGBw(2)/Lw);
    RGBp(:,3) = RGB(:,3)./(RGBw(3)/Lw);
else
    RGBwp(1) = RGBw(1)/(RGBw(2)/Lw);
    RGBwp(2) = RGBw(2)/(RGBw(2)/Lw);
    RGBwp(3) = RGBw(3)/(RGBw(2)/Lw);
    RGBp(:,1) = RGB(:,1)./(RGBw(2)/Lw);
    RGBp(:,2) = RGB(:,2)./(RGBw(2)/Lw);
    RGBp(:,3) = RGB(:,3)./(RGBw(2)/Lw);
end

LMSw =(MHM*RGBwp')';
LMS =(MHM*RGBp')';

rn = 0.57;
rsig = La;
LMSw(1) = LMSw(1)^rn/(LMSw(1)^rn + rsig^rn);
LMSw(2) = LMSw(2)^rn/(LMSw(2)^rn + rsig^rn);
LMSw(3) = LMSw(3)^rn/(LMSw(3)^rn + rsig^rn);

LMS(:,1) = LMS(:,1).^rn./(LMS(:,1).^rn + rsig.^rn);
LMS(:,2) = LMS(:,2).^rn./(LMS(:,2).^rn + rsig.^rn);
LMS(:,3) = LMS(:,3).^rn./(LMS(:,3).^rn + rsig.^rn);

crat = [40, 20, 1];
Aw = (crat(1)*LMSw(1) + crat(2)*LMSw(2) + crat(3)*LMSw(3))/sum(crat);
A = (crat(1).*LMS(:,1) + crat(2).*LMS(:,2) + crat(3).*LMS(:,3))./sum(crat);

A = A./Aw;

cn = 3.65;
csig = 0.65;
ca = 0.89;
cb = 0.24;
j = find(A<=cb);
k = find(A>=(ca+cb));
J = 100.*(-(A-cb).*(csig.^cn)./(A-cb-ca)).^(1/cn);
J(j) = 1; 
J(k) = 100;
J = J./100;
J = 100.*(mda.*(J-1)+1);
J(J<0) = 1;
J(J>100) = 100;
qa = 0.1308;
Q = J*(Lw.^qa);
arat = [11/11, -12/11, 1/11];
brat = [1/9, 1/9, -2/9];
a = arat(1).*LMS(:,1) + arat(2).*LMS(:,2) + arat(3).*LMS(:,3);
b = brat(1).*LMS(:,1) + brat(2).*LMS(:,2) + brat(3).*LMS(:,3);
Cscalar = 1358;
Cb = 0.6202;
Ca = 5.205;
C = Ca.*((Cscalar.*((a.^2+b.^2).^(1/2))).^Cb);
h=(180/pi)*atan2(real(b),real(a));
j=(h<0);h(j)=h(j)+360;
H = h;
j=find(h<=90 & h>=20.14); H(j)=(100.*(h(j)-20.14)./0.8)./(((h(j)-20.14)./0.8)+(90-h(j))/0.7);
j=find(h>=90 & h<=164.25); H(j)=(100+(100.*(h(j)-90)./0.7)./(((h(j)-90)./0.7)+(164.25-h(j))));
j=find(h>=164.25 & h<=237.53); H(j)=(200+(100.*(h(j)-164.25))./((h(j)-164.25)+((237.53-h(j))./1.2)));
j=find(h>=237.53 & h<=380.14); H(j)=(300+(100.*(h(j)-237.53)./1.2)./(((h(j)-237.53)./1.2)+(380.14-h(j))./0.8));
j=find(h<20.14); H(j)=(300+(100.*((h(j)+360)-237.53)./1.2)./((((h(j)+360)-237.53)./1.2)+(380.14-(h(j)+360))./0.8));
Ma = 0.11;
Mb = 0.61;
if (Lw==0) 
    Lw=1;
end
M = C.*(Ma*log10(Lw)+Mb);
s = 100.*(M./Q).^(1/2);
ac=C.*cos(h);
bc=C.*sin(h);
disp(sprintf('==================================='));
disp(sprintf('XYZw (input) : %.2f %.2f %.2f',XYZw(1), XYZw(2),XYZw(3) ));
disp(sprintf('La   (input) : %.2f cd/sqm (%3.2f%% )', La, 100*La/XYZw(1) ));
disp(sprintf('Media(input) : %s', media ));
disp(sprintf('==================================='));
end

% Inverse Model
function [XYZ] = xlrcami(J, C, h, Q, M, XYZw,La,media,CAT,sizeimg,mode)
    if strcmp(media,'hdr')
        mD = 1.0000;
    elseif strcmp(media,'trans')
        mD = 1.2175;
    elseif strcmp(media,'lcd')
        mD = 1.3374;
    elseif strcmp(media,'crt')
        mD = 1.4572;
    elseif strcmp(media,'paper')
        mD = 1.7526;
    end

    disp(sprintf('XYZw (output) : %.2f %.2f %.2f',XYZw(1), XYZw(2),XYZw(3) ));
    disp(sprintf('La   (output) : %.2f cd/sqm (%3.0f%% )', La, 100*La/XYZw(1) ));
    disp(sprintf('Media(output) : %s', media ));
    disp(sprintf('==================================='));

    % color transform matrices
    MCAT02=[0.7328,0.4296,-0.1624;
        -0.7036,1.6975,0.0061;
        0.0030,0.0136,0.9834];
    MCAT02i=[1.096124,-0.278869,0.182745;
        0.454369,0.473533,0.072098;
        -0.009628,-0.005698,1.015326];
    MHPE=[0.38971,0.68898,-0.07868;
        -.22981,1.1834,0.04641;
        0,0,1];
    MHPEi=[1.9102,-1.1121,0.2019;
        0.371,0.6291,-0;
        0,0,1;];
    MHM=MHPE*MCAT02i;
    MAabi =[1.0000, 0.3215, 0.2053;
            1.0000,-0.6351,-0.1860;
            1.0000,-0.1568,-4.4904;];

    Lw = XYZw(2);
    if (strcmp(mode,'QMh') || strcmp(mode,'QCh'))
        qa = 0.1308;
        J = Q/(Lw.^qa);
    end
    
    Jp = (1/mD).*(J./100-1)+1;
    RGBw =(MCAT02*XYZw')';
    if strcmp(CAT,'CAT')
        RGBwp(1) = RGBw(1)/(RGBw(1)/Lw);
        RGBwp(2) = RGBw(2)/(RGBw(2)/Lw);
        RGBwp(3) = RGBw(3)/(RGBw(3)/Lw);
    else
        RGBwp(1) = RGBw(1)/(RGBw(2)/Lw);
        RGBwp(2) = RGBw(2)/(RGBw(2)/Lw);
        RGBwp(3) = RGBw(3)/(RGBw(2)/Lw);
    end
    
    LMSw =(MHM*RGBwp')';
    
    rn = 0.57;
    LMSwp(1) = LMSw(1)^rn/(LMSw(1)^rn + La^rn);
    LMSwp(2) = LMSw(2)^rn/(LMSw(2)^rn + La^rn);
    LMSwp(3) = LMSw(3)^rn/(LMSw(3)^rn + La^rn);
    Aw = (40*LMSwp(1)+20*LMSwp(2)+LMSwp(3))/61;
    
    Aa = 0.89;
    Ab = 0.24;
    As = 0.65;
    An = 3.65;
    A = Aw.*(Aa.*Jp.^An./(Jp.^An+As^An)+Ab);
    
    if (strcmp(mode,'QMh') || strcmp(mode,'JMh'))
        Ma = 0.11;
        Mb = 0.61;
        C = M./(Ma*log10(Lw)+Mb);
    end
    
    Ca = 465.5;
    Cn = 0.62;
    CC = (C./Ca).^(1/Cn);
    r=(pi/180);
    a=cos(r*h).*CC;
    b=sin(r*h).*CC;
    LMSp = (MAabi*[A, a, b]')';
    Hn = 0.57;
    LMSp(LMSp<0)=0;
    LMS(:,1) = (-La^Hn.*LMSp(:,1)./(LMSp(:,1)-1)).^(1/Hn);
    LMS(:,2) = (-La^Hn.*LMSp(:,2)./(LMSp(:,2)-1)).^(1/Hn); 
    LMS(:,3) = (-La^Hn.*LMSp(:,3)./(LMSp(:,3)-1)).^(1/Hn);

    lmax = max(LMS(:,1));
    mmax = max(LMS(:,2));
    sclamp = max([lmax,mmax]);
    lms_s = LMS(:,3);
    lms_s(lms_s > sclamp) = sclamp;
    LMS(:,3) = lms_s;
    XYZp = (MHPEi*LMS')';
    XYZtmp = m2v(imresize(v2m(XYZp,sizeimg),0.3,'bicubic')); % ignoring noise, estimating white
    [CC, NN] = max(XYZtmp(:,2));
    XYZmax = XYZtmp(NN,:);
    nRGBw = RGBw./RGBw(2);
    Yscalar = XYZw(2)/XYZmax(2);
    if strcmp(CAT,'CAT')
        XYZtmp = m2v(imresize(v2m(XYZp,sizeimg),0.3,'bicubic')); % ignoring noise, estimating white
        [CC, NN] = max(XYZtmp(:,2));
        XYZmax = XYZtmp(NN,:);
        nRGBw = RGBw./RGBw(2);
        Yscalar = XYZw(2)/XYZmax(2);

        Mwb = Yscalar*[nRGBw(1),0,0;0,nRGBw(2),0;0,0,nRGBw(3)];
    else
        Yscalar = 1.0; % this part is newly added for the HDR video paper.
        Mwb = Yscalar*[1,0,0;0,1,0;0,0,1];
    end
    XYZ = (MCAT02i * Mwb * MCAT02 * XYZp')';
end