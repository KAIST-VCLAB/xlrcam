function [HDRi, fileinfo] = readhdr(filename)
    if nargin==0;
        [file_name,file_path] = uigetfile({'*.hdr','HDR'},'Choose an HDR image');
        source_images = strcat(file_path,file_name);   
        filename = char(source_images);
    end
    fid = fopen(filename,'r');

    tline = fgetl(fid);
    fileinfo.identifier = tline(3:end);

    tline = fgetl(fid);
    while ~isempty(tline)
        n = strfind(tline,'=');
        if ~isempty(n)
            vname = lower(tline(1:n(1)-1));
            vval = tline(n+1:end);
            fileinfo = setfield(fileinfo,vname,tline(n+1:end));
        end
        tline = fgetl(fid);
    end

    tline = fgetl(fid);
    fileinfo.Ysign = tline(1);
    [fileinfo.height,count,errmsg,nextindex] = sscanf(tline(4:end),'%d',1);
    fileinfo.Xsign = tline(nextindex+4);
    [fileinfo.width,count,errmsg,nextindex] = sscanf(tline(nextindex+7:end),'%d',1);

    HDRi = zeros(fileinfo.height, fileinfo.width, 3);
    [data, count] = fread(fid,inf,'uint8');
    fclose(fid);

    scanline_width = fileinfo.width;
    num_scanlines = fileinfo.height;

    if ((scanline_width < 8)||(scanline_width > 32767))
        HDRi = rgbe2float(reshape(data,fileinfo.width,fileinfo.height,4));
        return;
    end

    scanline_buffer = repmat(uint8(0),scanline_width,4);
    dp = 1;
    for scanline=1:num_scanlines
        if (data(dp) ~= 2) || (data(dp+1) ~= 2)
            rgbe(:,1) = data(1:4:end);
            rgbe(:,2) = data(2:4:end);
            rgbe(:,3) = data(3:4:end);
            rgbe(:,4) = data(4:4:end);

            float = rgbe2float(rgbe);
            img = reshape(float, fileinfo.width, fileinfo.height, 3);
            HDRi = transpose3(img);
            return;
        end

        if (bitshift(data(dp+2),8)+data(dp+3))~=scanline_width
            error('wrong scanline width');
        end
        dp = dp+4;
        for i=1:4 
            ptr = 1;
            while(ptr <= scanline_width)
                if (data(dp) > 128)
                    count = data(dp)-128;
                    if ((count == 0)||(count-1 > scanline_width - ptr)) 
                        warning('bad scanline data');
                    end
                    scanline_buffer(ptr:ptr+count-1,i) = data(dp+1);
                    dp = dp+2;
                    ptr = ptr+count;
                else % a non-run
                    count = data(dp);
                    dp = dp+1;
                    if ((count == 0)||(count-1 > scanline_width - ptr)) 
                        warning('bad scanline data');
                    end
                    scanline_buffer(ptr:ptr+count-1,i) = data(dp:dp+count-1);
                    ptr = ptr+count;
                    dp = dp+count;
                end
            end
        end
        HDRi(scanline,:,:) = rgbe2float(scanline_buffer);
    end 
end

function [rgb] = rgbe2float(rgbe)
    s = size(rgbe);
    rgbe = reshape(rgbe,prod(s)/4,4);
    rgb = zeros(prod(s)/4,3);
    l = find(rgbe(:,4)>0);
    rgb(l,:) = double(rgbe(l,1:3)).*repmat(2.^(double(rgbe(l,4))-128-8),1,3);
    rgb = reshape(rgb,[s(1:end-1),3]);
end
