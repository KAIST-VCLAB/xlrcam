% This is for siggraph asia Microscale video


num_images = 301;
start_frame_index = 0;
prefix_previous_job = 'leatherblue_';
scale = 0.01;

file_extension_HDR = '.exr';
out_directory = './';
dir_previous_job = 'exr';
outXYZw = d(65);
outLa = 10;
outmedia = 'trans';
CSpace = 'JMh';
CAT = 'noCAT';

addpath('/Users/minhkim/Dropbox/Research/codes/xlrcam/xlrcam_matlab');
save_TM_HDR = 1;
prefix_TM_HDR = 'leatherblue_';
file_extension = '.png';
dir_TM_HDR = 'png'

for i = 0:( num_images - 1 )
        frame_idx = start_frame_index + i;
        if frame_idx == 179
            continue
        end
        % Read images
        filename_to_read = sprintf( '%s%3.3d%s', prefix_previous_job, frame_idx, file_extension_HDR );
        fullfile_to_read = fullfile( out_directory, dir_previous_job, filename_to_read );
        %[ HDR_image, fileinfo] = readhdr( fullfile_to_read );
        HDR_image = scale.*readexr( fullfile_to_read );
        fprintf( '[HDR-TM] Processing %d th frame...( %d/%d )\n', frame_idx, i + 1, num_images );
        if i == 0 % we estimate adaptive parameters 
            [ HDR_TM_image, La, inXYZwhite, XYZ2ma ] = hlcamrep( HDR_image, outXYZw, outLa, outmedia, CSpace, CAT );
        else
            [ HDR_TM_image, La, inXYZwhite, XYZ2ma ] = hlcamrep( HDR_image, outXYZw, outLa, outmedia, CSpace, CAT, La, inXYZwhite, XYZ2ma );
        end
        if save_TM_HDR
            img_to_save = HDR_TM_image;
            file_name_to_save = sprintf( '%s%d%s', prefix_TM_HDR, frame_idx, file_extension );
            fullfile_to_save = fullfile( out_directory, dir_TM_HDR, file_name_to_save );
            imwrite( img_to_save, fullfile_to_save );
        end
    end