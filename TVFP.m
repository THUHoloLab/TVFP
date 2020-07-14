function [f,Out] = TVFP(b,samplingIndices,pupil,prms,f_true)

% ADM solves the ptychographic reconstruction based on total variation

% Inputs:
%  b        -- measured magnitude  
%  samplingIndices  -- sampling indices
%  pupil    -- pupil funciton
%  prms     -- parameters for algorithm
%  f_true	-- (optional) true image
%
% Outputs: 
%          
%  f      -- reconstructed image
%  Out    -- iteration information, e.g., iter number, relative errors,
%            function values, etc.
% -----------------------------------------------------------
% Jiachen Wu, created on Mar. 24, 2020 

%% set options
maxItr = prms.maxItr;
lambda = prms.lambda;
eta = prms.eta;
mu = prms.mu;
beta = prms.beta;

hROW = prms.hROW;
hCOL = prms.hCOL;

if ~isfield(prms,'tvtype')
    TVtype = 'isotropic';    % default TV type 
else
    TVtype = prms.tvtype;
end

switch TVtype
    case 'anisotropic'
        w_update = @(sh,sv) soft_anisotropic(sh,sv,eta/lambda);
    case 'isotropic'
        w_update = @(sh,sv) soft_isotropic(sh,sv,eta/lambda);
    otherwise
        error('Unrecognized TV type; it should be ''anisotropic'' or ''isotropic''');
end

b = sqrt(b);
[h,w,~] = size(b);

% set up the annoyomous helper functions
f_h = @(z) F_LENS2SENSOR(z(samplingIndices),pupil,h,w);
ft_h = @(z) F_SENSOR2LENS(z,samplingIndices,hROW,hCOL,h,w,pupil);

relchg_tol = prms.relchg_tol;
verbose = prms.verbose;       % screen display switch; turning it on slows the code down
prev = prms.prev;

%% initialization 
Hh = psf2otf([-1 1],[hROW,hCOL]);
Hh = fftshift(Hh);

HhT = conj(Hh);
Hv = psf2otf([-1;1],[hROW,hCOL]);
Hv = fftshift(Hv);

HvT = conj(Hv);

Ih = Hh.*HhT;
Iv = Hv.*HvT;

h = f_h(ones(hROW,hCOL));
sum_P = real((ft_h(h)));
denominator = mu*sum_P + eta*(Ih + Iv);

% Fb = fft2(b);
% Fb = fftshift(fftshift(Fb,1),2);
% x = ifft2(Fb);

f = imresize(b(:,:,round(end/2)),[hROW,hCOL]);
psi = fftshift(fft2(f));

% auxiliary variables
wh = 0;
wv = 0;
vh = 0;
vv = 0;
y = 0;
% P = double(P);

% initial SNR computation
if prev
   figure, im = imagesc(f);t = title('iteration = 0');
   axis image off tight
   colormap gray
   
%    figure, an = animatedline;
%    xlim([0,maxItr])
end

%% Main loop
for ii = 1:maxItr
        
    % ================================
    %  Begin Alternating Minimization

    % ----------------
    %   x-subprolem
    % ----------------
    fpsi = f_h(psi);
    s = fpsi - y/mu;
    x = sign(s).*(b + mu*abs(s))/(1 + mu); 

    % ----------------
    %   psi-subprolem
    % ----------------
%     y_inv = ifft2(y);
    numerator_1 = mu*ft_h(x + y/mu);            
    numerator_2 = eta*(HhT.*fft2(wh + vh/eta) + HvT.*fft2(wv + vv/eta));
    psi = (numerator_1 + numerator_2 + eps) ./ (denominator + eps);
    
    % ----------------
    %   w-subprolem
    % ----------------
    ipsi_h = ifft2(Hh.*psi);
    ipsi_v = ifft2(Hv.*psi);
    sh = ipsi_h - vh/eta;
    sv = ipsi_v - vv/eta;

    [wh,wv] = w_update(sh,sv);
                
    % --------------------
    % update dual variable
    % --------------------
    vh = vh + beta*eta*(wh - ipsi_h);
    vv = vv + beta*eta*(wv - ipsi_v);
    y = y + beta*mu*(x - fpsi);
   
    %  End Alternating Minimization
    % ================================

    % ----------------------------
    % check stopping criterion
    % ----------------------------
    f_prev = f;
    f = ifft2(fftshift(psi));
    
    relchg = norm(f-f_prev,'fro')/norm(f,'fro');
    
    if verbose && mod(ii,10) == 0
        fprintf('itr=%d relchg=%4.1e', ii, relchg);
        if exist('f_true','var')
            fprintf(' snr=%4.1f',snr(f_true,f - f_true)); 
        end
        fprintf('\n');
    end
    
    if prev
        set(im,'CData',abs(f));
        set(t,'String',['Iteration = ',num2str(ii)]);
        
%         regul = sqrt(abs(ipsi_h).^2 + abs(ipsi_v).^2);
%         regul = sum(regul,'all');
%         resid = abs(x - b).^2;
%         resid = sum(resid,'all');
%         objective = regul + lambda*resid/2;
%         addpoints(an,ii,double(resid));
        
        drawnow
    end

    if relchg < relchg_tol
        break;
    end
    
end

% outer
Out.iter = ii;
Out.relchg = relchg;
Out.psi = psi;

% if exist('f_true','var')
% Out.rmse = norm((f - f_true)/sqrt(hROW*hCOL),'fro');
% Out.SSIM = ssim(f,f_true);
% Out.PSNR = psnr(double(f),f_true);
% end

end

%% functions

function [wh,wv] = soft_anisotropic(sh,sv,eta)

    wh = sign(sh).*max(abs(sh) - 1/eta,0);
    wv = sign(sv).*max(abs(sv) - 1/eta,0);

end

function [wh,wv] = soft_isotropic(sh,sv,eta)

    s = sqrt(abs(sh).^2 + abs(sv).^2);
    wh = (sh./(s + eps)).*max(s - 1/eta,0);
    wv = (sv./(s + eps)).*max(s - 1/eta,0);

end