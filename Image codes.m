% read image and extract its components 
% read image
image = imread('peppers.png');
% red comp
red_component = image;
red_component(:,:,2)=0;
red_component(:,:,3)=0;
% green comp
green_component = image;
green_component(:,:,1)=0;
green_component(:,:,3)=0;
% blue comp 
blue_component = image;
blue_component(:,:,1)=0;
blue_component(:,:,2)=0;

%to display original image
figure;
imshow(image); title('original image');
%to display red component
figure;
imshow(red_component); title('red component'); 
 %to display green component
 figure;
imshow(green_component); title('green component');
%to display blue componen
 figure;
imshow(blue_component); title('blue component'); 

%%

% detecting the edges of the image 
image = imread('peppers.png'); 

red_comp = image(:, :, 1);
green_comp = image(:, :, 2);
blue_comp = image(:, :, 3);

% Convert to double 
red_comp = double(red_comp) / 255.0;
green_comp = double(green_comp) / 255.0;
blue_comp = double(blue_comp) / 255.0;

% Edge kernel 
edge_kernel_x = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
edge_kernel_y = [1, 2, 1; 0, 0, 0; -1, -2, -1];

% Convolution using conv2 

% red
edge_redx = conv2(red_comp, edge_kernel_x, 'same');
edge_redy = conv2(red_comp, edge_kernel_y, 'same');
edge_red = sqrt((edge_redx.^2)+(edge_redy.^2));

% green
edge_greenx = conv2(green_comp, edge_kernel_x, 'same');
edge_greeny = conv2(green_comp, edge_kernel_y, 'same');
edge_green = sqrt((edge_greenx.^2)+(edge_greeny.^2));

% blue
edge_bluex = conv2(blue_comp, edge_kernel_x, 'same');
edge_bluey = conv2(blue_comp, edge_kernel_y, 'same');
edge_blue = sqrt((edge_bluex.^2)+(edge_bluey.^2));

% combining
edge_image = cat(3, edge_red, edge_green, edge_blue);

% result
figure;
imshow(edge_image); title('Edge Image');

%%

% image sharping
image = imread('peppers.png'); 

red_comp = image(:, :, 1);
green_comp = image(:, :, 2);
blue_comp = image(:, :, 3);

% Convert to double 
red_comp = double(red_comp) / 255.0;
green_comp = double(green_comp) / 255.0;
blue_comp = double(blue_comp) / 255.0;

% sharping kernel 
sharping_kernel = [0, -1, 0; -1, 5, -1; 0, -1, 0];

% Convolution using conv2 

% red
sharp_red = conv2(red_comp, sharping_kernel, 'same');

% green
sharp_green = conv2(green_comp, sharping_kernel, 'same');

% blue
sharp_blue = conv2(blue_comp, sharping_kernel, 'same');

% combining
sharp_image = cat(3, sharp_red, sharp_green, sharp_blue);

% result
figure;
imshow(sharp_image); title('sharpened Image');

%%

% image bluring
image = imread('peppers.png');

red_comp = image(:, :, 1);
green_comp = image(:, :, 2);
blue_comp = image(:, :, 3);

% Convert to double
red_comp = double(red_comp) / 255.0;
green_comp = double(green_comp) / 255.0;
blue_comp = double(blue_comp) / 255.0;

% blurring kernel
% to make blur kernel we will take the average
bluring_kernel = [1,1,1,1,1 ; 1,1,1,1,1 ; 1,1,1,1,1 ; 1,1,1,1,1 ; 1,1,1,1,1] / 25 ;

% Convolution using conv2

% red
blur_red = conv2(red_comp, bluring_kernel, 'same');

% green
blur_green = conv2(green_comp, bluring_kernel, 'same');

% blue
blur_blue = conv2(blue_comp, bluring_kernel, 'same');

% combining
blur_image = cat(3, blur_red, blur_green, blur_blue);

% result
figure;
imshow(blur_image); title('Blurred Image');

%%

% image motion bluring
image = imread('peppers.png');

red_comp = image(:, :, 1);
green_comp = image(:, :, 2);
blue_comp = image(:, :, 3);

% Convert to double
red_comp = double(red_comp) / 255.0;
green_comp = double(green_comp) / 255.0;
blue_comp = double(blue_comp) / 255.0;

% Motion blurring kernel
Mbluring_kernel_lenght = 15; % Length of the motion blur
Mbluring_kernel = ones(1, Mbluring_kernel_lenght) / Mbluring_kernel_lenght;

% Convolution using conv2

% red
Mblur_red = conv2(red_comp, Mbluring_kernel, 'same');

% green
Mblur_green = conv2(green_comp, Mbluring_kernel, 'same');

% blue
Mblur_blue = conv2(blue_comp, Mbluring_kernel, 'same');

% combining
blur_image = cat(3, blur_red, blur_green, blur_blue);

% result
figure;
imshow(blur_image); title('Blurred Image');

%%

% restore original image from motion blurred one 

% firstly image motion bluring
image = imread('peppers.png');

red_comp = image(:, :, 1);
green_comp = image(:, :, 2);
blue_comp = image(:, :, 3);

% Convert to double
red_comp = double(red_comp) / 255.0;
green_comp = double(green_comp) / 255.0;
blue_comp = double(blue_comp) / 255.0;

% Motion blurring kernel
Mbluring_kernel_lenght = 15; % Length of the motion blur
Mbluring_kernel = ones(1, Mbluring_kernel_lenght) / Mbluring_kernel_lenght;

% Convolution using conv2

% red
Mblur_red = conv2(red_comp, Mbluring_kernel);

% green
Mblur_green = conv2(green_comp, Mbluring_kernel);

% blue
Mblur_blue = conv2(blue_comp, Mbluring_kernel);

% combining to get blurred image
Mblur_image = cat(3, Mblur_red, Mblur_green, Mblur_blue);


% now restoring thw original

% Compute the 2D Fourier transform of the motion blur kernel
Mblur_Kernel_fft = fft2(Mbluring_kernel, size(Mblur_image, 1) ,size(Mblur_image, 2));

% Compute the 2D Fourier transform of the blurred image
Mblur_image_fft = fft2(Mblur_image);

% Add a small constant to the denominator to avoid division by zero
epsilon = 0.0001 ;
Mblur_Kernel_fft = Mblur_Kernel_fft + epsilon;

% Perform deconvolution in the frequency domain
restored_Image_fft = Mblur_image_fft ./(Mblur_Kernel_fft + epsilon);

% Compute the inverse Fourier transform to get the restored image
restored_Image = abs(ifft2(restored_Image_fft));

% Crop the restored image to match the size of the original image
restored_Image = restored_Image(1:size(image, 1), 1:size(image, 2), :);

% Display the original and restored images
figure;
subplot(1,2,1); imshow(Mblur_image); title(' Motion Blurred Image');
subplot(1, 2, 2); imshow(restored_Image), title('Restored Image');