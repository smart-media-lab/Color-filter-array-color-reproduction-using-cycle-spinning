function CFA_IJEC_2010_main
%
% This is a demo program of color filter array (CFA) image color reproduction.
% J. Tian, W. Yu and L. Ma, "Color filter array color reproduction using 
% cycle-spinning," International Journal of Electronics and Communications,
% Vol. 64, No. 6, Jun. 2010, pp. 584-587.
%

close all; clear all; clc;
   
infilename = 'test.bmp';
wavetype = 'bior4.4'; %wavelet basis
nCS = 6; %maximum shift range of cycle spinning

img_orig = imread(infilename);
img_orig = im2double(img_orig);

% Generate the image in Bayer pattern
img_bayer = func_full2bayer(img_orig);

% Perform image demosaicking
img_out_interp_CFA_CS = zeros(size(img_orig));
for i = 0:2:nCS
    for j = 0:2:nCS
        img_temp = func_cyclespin2(img_bayer, i, j);
        out_temp = func_interp_CFA(img_temp, wavetype);

        dout = zeros(size(img_orig));
        for k =1:3
            dout(:,:,k) = func_cyclespin2(out_temp(:,:,k), size(out_temp(:,:,k),1)-i,size(out_temp(:,:,k),2)-j);
        end        
        img_out_interp_CFA_CS = img_out_interp_CFA_CS + dout;
    end
end

img_out_interp_CFA_CS = img_out_interp_CFA_CS./((nCS/2+1)^2);

% Write the output image
imwrite(img_out_interp_CFA_CS, 'test_rec.bmp', 'bmp');

% Calculate the PSNR performance
fprintf('PSNR (averaged over RGB channels) is %.2f dB.\n', func_psnr_color(img_orig.*255, img_out_interp_CFA_CS.*255));

%--------------------------------------------------------------------------
%----------------------------- Inner Function -----------------------------
%--------------------------------------------------------------------------
function result = func_full2bayer(RGB)
% generate Bayer data from the input image with full color components.
% G R
% B G

result = zeros(size(RGB,1),size(RGB,2));
result(1:2:end,2:2:end) = RGB(1:2:end,2:2:end,1); % Red
result(1:2:end,1:2:end) = RGB(1:2:end,1:2:end,2); % Green
result(2:2:end,2:2:end) = RGB(2:2:end,2:2:end,2); % Green
result(2:2:end,1:2:end) = RGB(2:2:end,1:2:end,3); % Blue

%--------------------------------------------------------------------------
%----------------------------- Inner Function -----------------------------
%--------------------------------------------------------------------------
function s = func_cyclespin2(x,i,j)
%generate the cycle spinning version of the input x according to its row
%and column shifts, i and j, respectively.
[l1,l2] = size(x);
z = x((l1+1-i):l1,:);
z(i+1:l1,:) = x(1:(l1-i),:);
s = z(:,(l2+1-j):l2);
s(:,j+1:l2) = z(:,1:(l2-j));

%--------------------------------------------------------------------------
%----------------------------- Inner Function -----------------------------
%--------------------------------------------------------------------------
function result = func_interp_CFA(X,wavtype)
%interpolate the CFA image
% G R
% B G

%%%%% Initial interpolation of the color channels
RGB = func_interp_bilinear(X);

%%%%% Initial Red, green, and blue channel assignments
R = RGB(:,:,1); 
G = RGB(:,:,2); 
B = RGB(:,:,3); 

%%%%%%%%%%%%%%%%%
R11=R(1:2:end,1:2:end);
R12=R(1:2:end,2:2:end); %%
R21=R(2:2:end,1:2:end);
R22=R(2:2:end,2:2:end);
G11=G(1:2:end,1:2:end); %%
G12=G(1:2:end,2:2:end);
G21=G(2:2:end,1:2:end);
G22=G(2:2:end,2:2:end); %%
B11=B(1:2:end,1:2:end);
B12=B(1:2:end,2:2:end);
B21=B(2:2:end,1:2:end); %%
B22=B(2:2:end,2:2:end);

%%subband decomposition and swap 
G12 = func_interp_CFA_subband_swap(G12,R12,wavtype);
B12 = func_interp_CFA_subband_swap(B12,R12,wavtype);
G21 = func_interp_CFA_subband_swap(G21,B21,wavtype);
R21 = func_interp_CFA_subband_swap(R21,B21,wavtype);
R11 = func_interp_CFA_subband_swap(R11,G11,wavtype);
B11 = func_interp_CFA_subband_swap(B11,G11,wavtype);
R22 = func_interp_CFA_subband_swap(R22,G22,wavtype);
B22 = func_interp_CFA_subband_swap(B22,G22,wavtype);

%%subband decomposition and swap 
R(1:2:end,1:2:end)=R11;
R(1:2:end,2:2:end)=R12; %%
R(2:2:end,1:2:end)=R21;
R(2:2:end,2:2:end)=R22;
G(1:2:end,1:2:end)=G11; %%
G(1:2:end,2:2:end)=G12;
G(2:2:end,1:2:end)=G21;
G(2:2:end,2:2:end)=G22; %%
B(1:2:end,1:2:end)=B11;
B(1:2:end,2:2:end)=B12;
B(2:2:end,1:2:end)=B21; %%
B(2:2:end,2:2:end)=B22;

%%%%% Output the image...
result(:,:,1) = R;
result(:,:,2) = G;
result(:,:,3) = B;
result(result<0) = 0;
result(result>1) = 1;

%--------------------------------------------------------------------------
%----------------------------- Inner Function -----------------------------
%--------------------------------------------------------------------------
function result = func_interp_CFA_subband_swap(A,B,wavtype)
%inner function used in func_interp_CFA

[AL,AH,AV,AD] = dwt2(A,wavtype);
[BL,BH,BV,BD] = dwt2(B,wavtype); 
BH=BH+(mean(AH(:))-mean(BH(:)));
BV=BV+(mean(AV(:))-mean(BV(:)));
BD=BD+(mean(AD(:))-mean(BD(:)));
result = idwt2(AL,BH,BV,BD,wavtype);

%--------------------------------------------------------------------------
%----------------------------- Inner Function -----------------------------
%--------------------------------------------------------------------------
function [RGB] = func_interp_bilinear(X)

% G R
% B G

% masks for Bayer pattern
Mr  = zeros(size(X)); Mr(1:2:end,2:2:end) = 1;   % 0:not red,  1:red
Mb  = zeros(size(X)); Mb(2:2:end,1:2:end) = 1;   % 0:not blue, 1:blue
Mg1 = zeros(size(X)); Mg1(1:2:end,1:2:end) = 1;  % 0: not G1,  1:G1
Mg2 = zeros(size(X)); Mg2(2:2:end,2:2:end) = 1;  % 0: not G2,  1:G2

% Green
Hg = [0 1/4 0; 1/4 1 1/4; 0 1/4 0];
G = conv2((Mg1+Mg2).*X,Hg,'same');

% Red/Blue
Hr = [1/4 1/2 1/4; 1/2 1 1/2; 1/4 1/2 1/4];
R = conv2(Mr.*X,Hr,'same');
B = conv2(Mb.*X,Hr,'same');

RGB(:,:,1) = R;
RGB(:,:,2) = G;
RGB(:,:,3) = B;

RGB(RGB<0) = 0;
RGB(RGB>1) = 1;

%--------------------------------------------------------------------------
%----------------------------- Inner Function -----------------------------
%--------------------------------------------------------------------------
function result = func_psnr_color(img1,img2)
%calculate the psnr value of two images, by averaging over R, G and B
%channels.
img1 = double(img1);
img2 = double(img2);
result = 0;
for i=1:3
    f = img1(:,:,i);
    g = img2(:,:,i);
    Q = 255; MSE = 0;
    [M,N] = size(f);
    h = f - g;
    MSE = sum(sum(h .* h));
    MSE=MSE/M/N;
    result=10*log10(Q*Q/MSE)+result;
end

result = result/3;