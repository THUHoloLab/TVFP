function [samplingIndices,pupil,hROW,hCOL] = getSampling(opts)

ROW = opts.imHeight;
COL = opts.imWidth;
% Tx = opts.nX;
% Ty = opts.nY;
% apShift = opts.apertureShift;
cc = opts.cc(:,1:2);
pupilType = opts.pupilType;
% samplingPattern = opts.samplingPattern;

Num = size(cc,1);

if ~isfield(opts,'upSampling')
    upSampling = true;    % default upSampling 
else
    upSampling = opts.upSampling;
end

% create the pupil
switch lower(pupilType)
    case 'circle'
        apDia = opts.apDia;
        pupil = imresize(padarray(fspecial('disk', apDia/2), floor((ROW - apDia)/2)*[1 1]), [ROW COL],'bilinear');
        pupil(pupil~=0)=1;
%         [x,y] = meshgrid((1:ROW)-ROW/2-1,(1:COL)-COL/2-1);
%         pupil = x.^2 + y.^2 < apDia^2/4;
%     case 'aberration'
%         % Zernike polynomial
%         [x,y] = meshgrid((1:ROW)-ROW/2-1,(1:COL)-COL/2-1);
%         a = opts.Zernike;
%         pupil = zeros(ROW,COL);
%         order = (sqrt(1 + 8*length(a)) - 1)/2 -1;
%         N = [];
%         M = [];
%         for ii = 0:order
%             N = [N ii*ones(1,ii+1)];
%             M = [M -ii:2:ii];
%         end
%         [theta,rho] = cart2pol(x,y);
%         rho = rho / (0.5*apDia);
%         is_in_circle = find(rho(:) <= 1);
%         Z = zernfun(N,M,rho(is_in_circle),theta(is_in_circle));
%         aberration = exp(1i*1*Z*a);
%         pupil(is_in_circle) = aberration;
    case 'square'
        pupil = ones([ROW COL]);
    case 'custom'
        apDia = opts.cc(:,3);
        lambda = opts.cc(:,4);
        pupil = zeros(ROW,COL,Num);
        [xc,yc] = meshgrid((1:ROW)-ROW/2-1,(1:COL)-COL/2-1);
        r = sqrt(xc.^2 + yc.^2);
        for ii = 1:Num
            pupil(:,:,ii) = 1./(1 + exp(lambda(ii)*(r-apDia(ii)/2)));
        end
    otherwise
        error('Pupil type not supported');
end

pupil = single(pupil);

% create the sampling indices
if upSampling
    span = max(cc,[],1) - min(cc,[],1);
    min_x = min(cc(:,1));
    min_y = min(cc(:,2));
    hROW = ROW+floor(span(2));
    hCOL = COL+floor(span(1));
    update_ind = @(yi,xi) shift(yi,xi,ROW,COL);
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
if any(cc - floor(cc),'all')
    subpixel = true;
    pupil_mask = zeros([ROW COL Num],'single');
else
    subpixel = false;
end

samplingIndices = zeros(ROW,COL,Num);
count = 0;
% for tr = 1:Ty
    for k = 1:Num
%         rind = repmat(VEC((1:ROW)+(tc-1)*apShift), ROW,1); % no partial shifts allowed
%         cind = VEC(repmat((1:COL)+(tr-1)*apShift+1, COL,1)); % no partial shifts allowed
%         if ~samplingPattern(tc,tr)
%             continue; % no data was captured here, skip
%         end
        count = count+1;
%         yi = (tc-1)*apShift;
%         xi = (tr-1)*apShift;
        yi = cc(k,2) - min_y;
        xi = cc(k,1) - min_x;
        
        [rind, cind] = update_ind(yi,xi);
        samplingIndices(:,:,count) = reshape(sub2ind([hROW,hCOL],rind,cind),ROW,COL);
        
        if subpixel
            aper = imtranslate(pupil,[yi-floor(yi) xi-floor(xi)],'linear');
            if size(aper,3) == 1 
            pupil_mask(:,:,count) = aper;
            else
            pupil_mask(:,:,count) = aper(:,:,count);
            end
        end
    end
% end

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