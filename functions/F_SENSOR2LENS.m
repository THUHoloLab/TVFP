function x = F_SENSOR2LENS(y,samplingIndices,hROW,hCOL,ROW,COL,pupil)
suby = VEC(samplingIndices);

FY = fft2(y) .* conj(pupil);
subp = find(FY);
val = [FY(subp);0];
subs = [suby(subp);hROW * hCOL];

x = reshape(accumarray(subs,val),hROW,hCOL);
% x = fftshift(x);
end