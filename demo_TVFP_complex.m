clearvars -except PSNR RMSE; 
clc; close all

addpath('./functions')

im = im2double(imread('./images/house256.png'));
phi = im2double(imread('./images/peppers256.png'));
im = im.*exp(1i*(0.5*pi*phi-pi/4));     % phase range [-pi/4,pi/4]

[h,w] = size(im);
N = 9;
apDia = 40;
spacing = 25;
r = apDia/2;
d = apDia - spacing;
p = 2*acos(spacing/apDia)/pi - d*sqrt(r^2 - (d/2)^2)/(pi*r^2);
samplingPattern = ones(N);

ROW = 64;           % LR image size
COL = 64;           % LR image size

opts = struct();
opts.imHeight = ROW;
opts.imWidth = COL;
opts.hROW = h;      % HR image size
opts.hCOL = w;      % HR image size
opts.nX = N;
opts.nY = N;
opts.apertureShift = spacing;
opts.apDia = apDia;
opts.pupilType = 'circle';  % 'circle'or'aberration'
opts.samplingPattern = samplingPattern;
% opts.upSampling = false;

% foward model
[samplingIndices,pupil,hROW,hCOL] = getSampling(opts);

X = fftshift(fft2(im)); % need to center the FFT
X = padarray(X,floor([(hROW-h)/2 (hCOL-w)/2]));
y = F_LENS2SENSOR(X(samplingIndices),pupil,ROW,COL);
I = abs(y).^2;

% add noise
SNR = inf;
I = addNoise(I,SNR);
I(I<0) = 0; % input cannot be negative (avoid noise causing a negative signal)

% figure,imshow(log1p(imtile(abs(y))),[]);

% parameters for TVFP
prms.maxItr = 1000;
prms.hROW = hROW;
prms.hCOL = hCOL;
prms.ROW = ROW;
prms.COL = COL;
prms.lambda = 20e-03;
prms.mu = 0.6;% high value for the low SNR data  
prms.eta = 10;
prms.beta = 0.6;

prms.relchg_tol = 1e-05;
prms.verbose = 0;           % display the relative change
prms.prev = 1;              % preview the iterative image

tic
[x,out] = TVFP(I,samplingIndices,pupil,prms,im);
toc

LRImage = imresize(sqrt(I(:,:,ceil(end/2))),[hROW,hCOL]);
figure('position',[100,100,1024,384]),
subplot(131),imshow(abs(LRImage),[]);title('LR image')
subplot(132),imshow(abs(x),[]);title('Recovered amplitude')
subplot(133),imshow(angle(-x),[]);title('Recovered phase')

% %---evaluation
% % x = exp(-1i*pi/4)*x/sign(x(1));   % phase calibration
% Fx = fftshift(fft2(fftshift(x)));
% Fx = Fx((hROW-h)/2+1:end - (hROW-h)/2,(hCOL-w)/2+1:end - (hCOL-w)/2);
% x = fftshift(ifft2(fftshift(Fx)));
% M1 = norm((x - im)/sqrt(hROW*hCOL),'fro');
% M2 = norm((-x - im)/sqrt(hROW*hCOL),'fro');
% M = min(M1,M2);
% psn = 20*log10(1./M);
% PSNR = [PSNR psn];