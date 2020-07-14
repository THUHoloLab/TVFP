clearvars; 
clc; close all

addpath('./functions')

im = im2double(imread('./images/house256.png'));
phi = im2double(imread('./images/peppers256.png'));
im = im.*exp(1i*(pi*phi/2 - pi/4));

[h,w] = size(im);
N = 5;
apDia = 40;
spacing = 30;
samplingPattern = ones(N);
ROW = h - floor((N - 1)*spacing);
COL = w - floor((N - 1)*spacing);
% ROW = h;
% COL = w;

opts = struct();
opts.imHeight = ROW;
opts.imWidth = COL;
opts.nX = N;
opts.nY = N;
opts.apertureShift = spacing;
opts.apDia = apDia;
opts.pupilType = 'aberration';  % 'circle','aberration'
opts.Zernike = [ 0.03;-0.44; 0.20;-0.21;-0.52;...
                -0.26;-0.14;-0.08;-0.13; 0.05;...
                -0.23; 0.17;-0.38; 0.18; 0.06];     % Zernike coefficient
opts.samplingPattern = samplingPattern;
% opts.upSampling = false;

% foward model
[samplingIndices,pupil,hROW,hCOL] = getSampling(opts);
X = fftshift(fft2(im)); % need to center the FFT
% X = padarray(X,floor([(hROW-ROW)/2 (hCOL-COL)/2]));
y = F_LENS2SENSOR(X(samplingIndices),pupil,ROW,COL);
I = abs(y).^2;

% figure,imshow(log1p(imtile(abs(y))),[]);

% parameters for eTVFP
prms.maxItr = 500;
prms.hROW = hROW;
prms.hCOL = hCOL;
prms.ROW = ROW;
prms.COL = COL;

prms.lambda1 = 1e-04;	%0.0001;  
prms.mu1 = 0.001;       %0.001;    
prms.eta = 0.01;        %0.01;    

prms.lambda2 = 8;       %8; 
prms.mu2 = 0.0007;      %0.0007;
prms.xi = 10;           %10;
prms.gamma = 500;       %500;

prms.beta = 0.1;
% prms.tvtype = 'anisotropic';

prms.relchg_tol = 1e-05;
prms.verbose = 0;       % display the relative change
prms.prev = 1;          % preview the iterative image

% se = strel('disk',32);
% sup = imdilate(abs(pupil),se);
sup = abs(pupil);
% sup = ones(ROW,COL);

tic
[x,P,out] = eTVFP(I,samplingIndices,sup,prms,im,pupil);
toc

LRImage = imresize(sqrt(I(:,:,ceil(end/2))),[hROW,hCOL]);
figure('position',[100,100,1024,384]),
subplot(131),imshow(abs(LRImage),[]);title('LR image')
subplot(132),imshow(abs(x),[]);title('Recovered amplitude')
subplot(133),imshow(angle(x),[]);title('Recovered phase')
figure,imagesc(angle(pupil(ROW/2-31:ROW/2+32,COL/2-31:COL/2+32)),[-1 1])
axis equal off
colormap(jet)
title('True pupil phase')
figure,imagesc(angle(P(ROW/2-31:ROW/2+32,COL/2-31:COL/2+32)),[-1 1])
axis equal off
colormap(jet)
title('Recovered pupil phase')