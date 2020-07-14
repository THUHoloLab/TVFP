clearvars; 
clc; close all

addpath('./functions')
load('./realData/USAF_7x7.mat');

apDia = 63;

% get the Fourier domain sampling pattern
[h,w,N] = size(ims);

% samplingPattern = ones(N,1);
samplingPattern = zeros(N,1);
% samplingPattern([4 11 17 18 19 22 23 24 25 ...
%     26 27 28 31 32 33 39 46]) = 1;
samplingPattern([4 18 22 24 25 26 28 32 46]) = 1;
% samplingPattern([17 18 19 24 25 26 31 32 33]) = 1;
% ims = reshape(ims,h,w,[]);

ind = find(~samplingPattern);
ims(:,:,ind) = [];
cc(ind,:) = [];

% get the sampling indices and pupil
opts.imHeight = h; opts.imWidth = w;
opts.nX = 7; opts.nY = 7;
opts.samplingPattern = samplingPattern;
opts.cc = cc;
opts.apertureShift = spacing; 
opts.apDia = apDia;
opts.pupilType = 'circle';%'circle' 'custom'
[samplingIndices,pupil,hROW,hCOL] = getSampling(opts);

% parameters for TVFP
prms.maxItr = 500;
prms.ROW = h;
prms.COL = w;
prms.hROW = hROW;
prms.hCOL = hCOL;
prms.lambda1 = 10;%8;
prms.mu1 = 0.8;%1;     
prms.eta = 4;%4;     % relaxing

prms.lambda2 = 20;%20; 
prms.mu2 = 0.0007;%0.0007;
prms.xi = 100;%100;
prms.gamma = 1;%1;

prms.beta = 0.05;%0.05
% prms.tvtype = 'anisotropic';

prms.relchg_tol = 1e-05;
prms.verbose = 0;
prms.prev = 1;

% [x,y] = meshgrid(-150:149);
% r = sqrt(x.^2 + y.^2);
% pupil = 1./(1 + exp(0.5*(r-apDia/2)));

tic
[x,out] = eTVFP_multi_pupil(ims,samplingIndices,pupil,prms);
toc

% show the results
centerView = imresize(ims(:,:,ceil(end/2)),[hROW,hCOL],'bilinear');
figure,imshowpair(flipud(sqrt(centerView)),flipud(abs(x)),'montage');
