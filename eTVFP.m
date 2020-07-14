function [f,P,Out] = eTVFP(b,samplingIndices,sup,prms,f_true,pupil)

% ADM solves the ptychographic reconstruction 
% with pupil recovery based on total variation

% Inputs:
%  b        -- measured magnitude  
%  samplingIndices  -- sampling indices
%  sup    -- support of pupil funciton
%  prms     -- parameters for algorithm
%  f_true	-- (optional) true image
%  pupil	-- (optional) true pupil funciton
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
lambda1 = prms.lambda1;
lambda2 = prms.lambda2;
eta = prms.eta;
mu1 = prms.mu1;
mu2 = prms.mu2;
xi = prms.xi;
gamma = prms.gamma;
beta = prms.beta;

hROW = prms.hROW;
hCOL = prms.hCOL;
ROW = prms.ROW;
COL = prms.COL;

if ~isfield(prms,'tvtype')
    TVtype = 'isotropic';    % default TV type 
else
    TVtype = prms.tvtype;
end

switch TVtype
    case 'anisotropic'
        w_update = @(sh,sv) soft_anisotropic(sh,sv,eta/lambda1);
    case 'isotropic'
        w_update = @(sh,sv) soft_isotropic(sh,sv,eta/lambda1);
    otherwise
        error('Unrecognized TV type; it should be ''anisotropic'' or ''isotropic''');
end

b = sqrt(b);
[h,w,~] = size(b);

% set up the annoyomous helper functions
f_h = @(z,P) F_LENS2SENSOR(z(samplingIndices),P,h,w);
ft_h = @(z,P) F_SENSOR2LENS(z,samplingIndices,hROW,hCOL,h,w,P);

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

Kh = psf2otf([-1 1],[ROW,COL]);
% Kh = fftshift(Kh);
KhT = conj(Kh);
Ph = Kh.*KhT;
% Ph = fftshift(Ph);
Kv = psf2otf([-1;1],[ROW,COL]);
KvT = conj(Kv);
Pv = Kv.*KvT;

denominator = eta*(Ih + Iv);

% Fb = fft2(b);
% Fb = fftshift(fftshift(Fb,1),2);
% x = ifft2(Fb);

f = imresize(b(:,:,round(end/2)),[hROW,hCOL]);
f = f/max(f(:));
psi = fftshift(fft2(f));
P = sup;

% auxiliary variables
wh = 0;
wv = 0;
mh = 0;
mv = 0;
vh = 0;
vv = 0;
uh = 0;
uv = 0;
y = 0;
z = 0;

rmse_f = [];
rmse_P = [];

% initial plot window
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
    fpsi = f_h(psi,P);
    s = fpsi - y/mu1;
    x = sign(s).*(b + mu1*abs(s))/(1 + mu1); 
    
    % ----------------
    %   psi-subprolem
    % ----------------
%     y_inv = ifft2(y);
    numerator_1 = mu1*ft_h(x + y/mu1,P);            
    numerator_2 = eta*(HhT.*fft2(wh + vh/eta) + HvT.*fft2(wv + vv/eta));
    h = f_h(ones(hROW,hCOL),P);
    sum_P = ft_h(h,P);
    psi = (numerator_1 + numerator_2 + eps) ./ (mu1*sum_P + denominator + eps);
    
    % ----------------
    %   P-subprolem
    % ----------------
    Spsi = psi(samplingIndices);
    Fh_num = mu2*sum(conj(Spsi).*fft2(x + y/mu2),3);
    Fh_denom = mu2*sum(abs(Spsi).^2,3) + gamma*ones(ROW,COL);
    for kk = 1:5
        
        Fh = (Fh_num + gamma*P + z)./Fh_denom;
        FP = (xi*(KhT.*fft2(mh + uh/xi) + KvT.*fft2(mv + uv/xi)) + gamma*fft2(Fh - z/gamma))./...
             (xi*(Ph + Pv) + gamma*ones(ROW,COL));
        P = ifft2(FP);
        P(abs(P) > 1) = sign(P(abs(P) > 1));
        P = P.*sup;
        
        % m-subprolem
        DP_h = ifft2(Kh.*fft2(P));
        DP_v = ifft2(Kv.*fft2(P));
        th = DP_h - uh/xi;
        tv = DP_v - uv/xi;

        [mh,mv] = soft_isotropic(th,tv,xi/lambda2);
        uh = uh + beta*xi*(mh -  DP_h);
        uv = uv + beta*xi*(mv -  DP_v);

        z = z + beta*gamma*(P - Fh);
        
    end

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
    y = y + beta*mu1*(x - fpsi);
   
    %  End Alternating Minimization
    % ================================

    % ----------------------------
    % check stopping criterion
    % ----------------------------
    f_prev = f;
    f = ifft2(fftshift(psi));
    
    relchg = norm(f-f_prev,'fro')/norm(f,'fro');
    rmse_f = [rmse_f, norm((f - f_true)/sqrt(hROW*hCOL),'fro')];
    rmse_P = [rmse_P, norm((P - pupil)/sqrt(ROW*COL),'fro')];
    
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
%         objective = regul + lambda1*resid/2;
%         addpoints(an,ii,double(resid));
        
        drawnow
    end

    if relchg < relchg_tol
        break;
    end
    
end

% outer
% rmse = sqrt(norm(f - f_true,'fro')/(m*n));
Out.iter = ii;
Out.relchg = relchg;
Out.rmse_f = rmse_f;
Out.rmse_P = rmse_P;
% Out.SSIM = ssim(f,f_true);
% Out.PSNR = psnr(f,f_true);
Out.psi = psi;

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