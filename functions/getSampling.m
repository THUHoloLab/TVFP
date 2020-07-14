function [samplingIndices,pupil,hROW,hCOL] = getSampling(opts)

ROW = opts.imHeight;
COL = opts.imWidth;
Tx = opts.nX;
Ty = opts.nY;
apShift = opts.apertureShift;
apDia = opts.apDia;
pupilType = opts.pupilType;
samplingPattern = opts.samplingPattern;

if ~isfield(opts, 'hROW')
    hROW_input  = -1;  
else
    hROW_input = opts.hROW;
end

if ~isfield(opts, 'hCOL')
    hCOL_input  = -1;  
else
    hCOL_input = opts.hCOL;
end

if ~isfield(opts,'upSampling')
    upSampling = true;    % default upSampling 
else
    upSampling = opts.upSampling;
end

% create the pupil
switch lower(pupilType)
    case 'circle'
        pupil = imresize(padarray(fspecial('disk', apDia/2), floor((ROW - apDia)/2)*[1 1]), [ROW COL],'bilinear');
        pupil(pupil~=0)=1;
%         [x,y] = meshgrid((1:ROW)-ROW/2-1,(1:COL)-COL/2-1);
%         pupil = x.^2 + y.^2 < apDia^2/4;
    case 'aberration'
        % Zernike polynomial
        [x,y] = meshgrid((1:ROW)-ROW/2-1,(1:COL)-COL/2-1);
        a = opts.Zernike;
        pupil = zeros(ROW,COL);
        order = (sqrt(1 + 8*length(a)) - 1)/2 -1;
        N = [];
        M = [];
        for ii = 0:order
            N = [N ii*ones(1,ii+1)];
            M = [M -ii:2:ii];
        end
        [theta,rho] = cart2pol(x,y);
        rho = rho / (0.5*apDia);
        is_in_circle = find(rho(:) <= 1);
        Z = zernfun(N,M,rho(is_in_circle),theta(is_in_circle));
        aberration = exp(1i*1*Z*a);
        pupil(is_in_circle) = aberration;
    case 'square'
        pupil = ones([ROW COL]);
    otherwise
        error('Pupil type not supported');
end
pupil = single(pupil);

% create the sampling indices
if upSampling
    hROW_cal = ROW+floor(apShift*(Ty-1));
    hCOL_cal = COL+floor(apShift*(Tx-1));
    hROW = max(hROW_cal,hROW_input);
    hCOL = max(hCOL_cal,hCOL_input);
    update_ind = @(yi,xi) shift(yi+(hROW-hROW_cal)/2,xi+(hCOL-hCOL_cal)/2,ROW,COL);
else
    hROW = ROW;
    hCOL = COL;
    dx = floor(apShift*(Tx-1)/2);
    dy = floor(apShift*(Ty-1)/2);
    update_ind = @(yi,xi) cyclic_shift(yi,xi,ROW,COL,dy,dx);
end

% % check to make sure hROW and hCOL are correctly sized
% if (apShift-floor(apShift)) > eps
%     if mod(hROW,2)==1
%         hROW = hROW+1;
%     end
%     if mod(hCOL,2)==1
%         hCOL = hCOL+1;
%     end
% end

% check to see if the shift is fractional
if (apShift - floor(apShift)) > eps
    subpixel = true;
    pupil_mask = zeros([size(pupil) sum(samplingPattern(:))],'single');
else
    subpixel = false;
end

samplingIndices = zeros(ROW,COL,sum(samplingPattern(:)));
count = 0;
for tr = 1:Ty
    for tc = 1:Tx
%         rind = repmat(VEC((1:ROW)+(tc-1)*apShift), ROW,1); % no partial shifts allowed
%         cind = VEC(repmat((1:COL)+(tr-1)*apShift+1, COL,1)); % no partial shifts allowed
        if ~samplingPattern(tc,tr)
            continue; % no data was captured here, skip
        end
        count = count+1;
        yi = (tc-1)*apShift;
        xi = (tr-1)*apShift;
%         rind = repmat(VEC((1:ROW)+floor(yi)), COL,1);
%         cind = VEC(repmat((1:COL)+floor(xi), ROW,1));
        [rind, cind] = update_ind(yi,xi);
        samplingIndices(:,:,count) = reshape(sub2ind([hROW,hCOL],rind,cind),ROW,COL);
        
        if subpixel
            aper = imtranslate(pupil,[yi-floor(yi) xi-floor(xi)],'linear');
            pupil_mask(:,:,count) = aper;
        end
    end
end

if subpixel
    pupil = pupil_mask;
end

end

%% functions
function [rind, cind] = shift(yi,xi,ROW,COL)
	rind = repmat(VEC((1:ROW)+floor(yi)), COL,1);
	cind = VEC(repmat((1:COL)+floor(xi), ROW,1));
end

function [rind, cind] = cyclic_shift(yi,xi,ROW,COL,dx,dy)
    rind = repmat(VEC(circshift(1:ROW,dy-floor(yi))), COL,1);
    cind = VEC(repmat(circshift(1:COL,dx-floor(xi)), ROW,1));
end