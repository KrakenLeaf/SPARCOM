function [ Flag ] = DisplaySPE( InFile )
% DISPLAYSPE: Display an SPE image file.
%
% Syntax:
% -------
% [ Flag ] = DisplaySPE( InFile )
%
% Input:
% ------
% InFile   -   name of the input *.spe file (include the suffix .spe).
%
% Output:
% -------
% OutFile  -   Returns 0 if successful, otherwise 1.
%

try
    % read input image (WinView .spe format)
    fid = fopen(InFile, 'r', 'l');
    header = fread(fid,2050,'*uint16');
    
    % parsing .spe header
    x_dim = double(header(22));
    y_dim = double(header(329));
    z_dim = double(header(724));
    
    %determining the data format in the source image
    img_mode = header(55);
    if img_mode == 0
        imgs = fread(fid, inf, 'float32');
    elseif img_mode == 1
        imgs = fread(fid, inf, 'uint32');
    elseif img_mode == 2
        imgs = fread(fid, inf, 'int16');
    elseif img_mode == 3
        imgs = fread(fid, inf, 'uint16');
    end
    
    fclose(fid);
    
    imgs = reshape(imgs, [x_dim y_dim z_dim]);
    %     imgs = (imgs - ccd_base_line) * photon_per_count;
    
    colormap('gray');
    imagesc(imgs);
    title(['Image: ' num2str(x_dim) ' X ' num2str(y_dim) ' pixels']);
    
    % Disable ticks
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    
    Flag = 0;
catch
    Flag = 1;
end



