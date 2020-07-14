clearvars
clc; close all

addpath('./functions')

im = im2double(imread('./images/resChart.tif'));

[h,w] = size(im);
N = 3;
apDia = 80;
spacing = 90;
samplingPattern = ones(N);

ROW = 128;              % LR image size
COL = 128;              % LR image size

opts = struct();
opts.imHeight = ROW;    
opts.imWidth = COL;
opts.hROW = h;          % HR image size
opts.hCOL = w;          % HR image size
opts.nX = N;
opts.nY = N;
opts.apertureShift = spacing;
opts.apDia = apDia;
opts.pupilType = 'circle';
opts.samplingPattern = samplingPattern;
% opts.upSampling = false;

% foward model
[samplingIndices,pupil,hROW,hCOL] = getSampling(opts);

X = fftshift(fft2(im)); % need to center the FFT
X = padarray(X,floor([(hROW-h)/2 (hCOL-w)/2]));
y = F_LENS2SENSOR(X(samplingIndices),pupil,ROW,COL);
I = abs(y).^2;

% parameters for TVFP
prms.maxItr = 1000;
prms.hROW = hROW;
prms.hCOL = hCOL;
prms.lambda = 0.001;
prms.mu = 0.02;    
prms.eta = 0.01;
prms.beta = 0.2;
% prms.tvtype = 'anisotropic';

prms.relchg_tol = 1e-05;
prms.verbose = 1;       % display the relative change
prms.prev = 1;          % preview the iterative image

tic
[x,out] = TVFP(I,samplingIndices,pupil,prms);
toc

% show the results
centerView = imresize(I(:,:,ceil(end/2)),[hROW,hCOL],'bilinear');
figure,imshowpair(sqrt(centerView),abs(x),'montage');
title('Left:Center input image   Right:Recovered phase')