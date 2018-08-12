function [ ImageOut ] = ImagePostProcessing( ImageIn, method, Params )
%IMAGEPOSTPROCESSING - perform image post processing for the super-resolution
%image. This function performs only visualization tricks.

%% Initialization
try
    ColorMap = Params.ColorMap;
catch
    ColorMap = 'hot';   
end
[M, N] = size(ImageIn);

%% Choose method
switch lower(method)
    case 'same'
        % No post rendering on the input image
        % ------------------------------------
        ImageOut = ImageIn;
    case 'sharp'
        % each non-zero pixel in the input image is given a value of 1
        % ------------------------------------------------------------
        Inds = find(ImageIn ~= 0);
        ImageIn(Inds) = 1;
        ImageOut = ImageIn;
    case 'gaussian'
        % Blur the input image with a gaussian
        % ------------------------------------
        % Take parameters
        try
            A = Params.A;
        catch
            A = 1;
        end
        try
            Sigma = Params.Sigma;
        catch 
            Sigma = 0.005;
        end 
        
        % Create a mesh for the gaussian
        X = linspace(-1, 1, M);
        Y = linspace(-1, 1, N);
        [X, Y] = meshgrid(X, Y);
        
        % Create the gaussian
        G = A*exp(-0.5*(X).^2/Sigma^2).*exp(-0.5*(Y).^2/Sigma^2);
        
        % Convolve image with the gaussian
        ImageOut = conv2(ImageIn,G,'same');
    otherwise
        % Unknown rendering method
        % ------------------------
        error('ImagePostProcessing: No valid rendering method.');
end

%% Visualize 
if (Params.ShowImage)
    figure;
    eval(['colormap ' ColorMap ';']);
    imagesc(ImageOut);
    title(['Super resolution image using method: ' method]);
end






















