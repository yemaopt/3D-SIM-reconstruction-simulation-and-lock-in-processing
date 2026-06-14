% Ye Ma, 2026.06.13
% Department of Biomedical Engineering, Johns Hopkins University

[ImgSIM, info] = A6_3DSIM_RawDatasetGeneration_StrongDegrade_func();

function [ImgSIM, info] = A6_3DSIM_RawDatasetGeneration_StrongDegrade_func(varargin)
% A6_3DSIM_RawDatasetGeneration_StrongDegrade
% A4-style 3D-SIM raw-data simulator with STRONG objective-fixed
% depth-dependent modulation washout. This keeps the A4's illumination pattern, scanning
% geometry, sample embedding, and PSF convolution, but replaces the ideal
% pattern P with:
%
%   P_eff = mean(P) + alpha(z) * (P - mean(P))
%
% so the local mean illumination is preserved while modulation contrast
% decays away from the focal plane.
%
% Optional phase-independent smooth background can also be added for stress
% testing Lock-in/background rejection.
%
% Example:
%   [ImgSIM,info] = A6_3DSIM_RawDatasetGeneration_StrongDegrade( ...
%       'zWashout',0.5,'alphaFloor',0, ...
%       'bgStrength',0.4,'saveFile','Img3DSIM_StrongDegrade.mat');

p = inputParser;
p.addParameter('sampleFile','Sample1000Tubule_3D_256Pixels.mat',@ischar);
p.addParameter('psfFile','PSF_256Pixels.mat',@ischar);
p.addParameter('saveFile','Img3DSIM_StrongDegrade.mat',@ischar);
p.addParameter('PixelSize',0.08,@isnumeric);          % um
p.addParameter('IncidentAngle',pi/3,@isnumeric);
p.addParameter('emWavelength',0.488,@isnumeric);     % um
p.addParameter('zWashout',0.5,@isnumeric);           % um
p.addParameter('alphaFloor',0.0,@isnumeric);        % residual modulation far from focus
p.addParameter('zScanRange',-64:63,@isnumeric);      % same as original A4
p.addParameter('bgStrength',0.4,@isnumeric);         % phase-independent haze strength
p.addParameter('bgSigma',8,@isnumeric);              % smoothing for haze object, pixels
p.addParameter('useSingle',true,@islogical);
p.addParameter('showFigures',true,@islogical);
p.parse(varargin{:});
par = p.Results;

S = load(par.sampleFile);
if ~isfield(S,'Sample'), error('sampleFile must contain variable Sample'); end
Sample = S.Sample;
P = load(par.psfFile);
if ~isfield(P,'PSF'), error('psfFile must contain variable PSF'); end
PSF = P.PSF;

L = size(Sample,1);
if any(size(Sample) ~= [L L L]) || any(size(PSF) ~= [L L L])
    error('Sample and PSF must both be cubic volumes of the same size.');
end

if par.useSingle
    Sample = single(Sample);
    PSF = single(PSF);
else
    Sample = double(Sample);
    PSF = double(PSF);
end

% Match A4: embed central half of the sample in zero volume to prevent
% circshift wrap-around from contaminating the scanned focal plane.
temp = zeros(L,L,L,'like',Sample);
slabStart = L/4 + 1;
slabEnd   = 3*L/4;
temp(:,:,slabStart:slabEnd) = Sample(:,:,slabStart:slabEnd);

coords = (-L/2):(L/2-1);
[xxx,yyy,zzz] = meshgrid(coords,coords,coords);
xxx = cast(xxx .* par.PixelSize,'like',Sample);
yyy = cast(yyy .* par.PixelSize,'like',Sample);
zzz = cast(zzz .* par.PixelSize,'like',Sample);

k  = 2*pi/par.emWavelength;
kr = k*sin(par.IncidentAngle);
kz = k*cos(par.IncidentAngle);
A1 = 1; A2 = 1; A3 = 1;
thetaSet = [0, pi/3, 2*pi/3];
phiVec = 2*pi/5.*(0:4);

alphaZ = par.alphaFloor + (1-par.alphaFloor).*exp(-(zzz./par.zWashout).^2);
alphaZ = cast(alphaZ,'like',Sample);

% Precompute PSF FFT once; original A4 recomputed this in the loop.
fPSF = fftshift(fftn(ifftshift(PSF)));

% Precompute all 15 effective patterns.
PatternEff = zeros(L,L,L,15,'like',Sample);
PatternIdealOneAngle = zeros(L,L,L,5,'like',Sample);
ii = 0;
for theta = thetaSet
    kx = kr*cos(theta);
    ky = kr*sin(theta);
    for pp = 1:5
        phi1 = phiVec(pp);
        phi2 = phi1/2;
        ii = ii + 1;
        E1 = A1.*exp(1j.*(kx.*xxx + ky.*yyy + phi1) + 1j.*kz.*zzz);
        E2 = A2.*exp(1j.*(k.*zzz + phi2));
        E3 = A3.*exp(-1j.*(kx.*xxx + ky.*yyy) + 1j.*kz.*zzz);
        Pideal = cast(real(abs(E1 + E2 + E3).^2),'like',Sample);
        Pmean = mean(Pideal(:));
        PatternEff(:,:,:,ii) = Pmean + alphaZ .* (Pideal - Pmean);
        if theta == thetaSet(1)
            PatternIdealOneAngle(:,:,:,pp) = Pideal;
        end
    end
end

ImgSIM = zeros(L,L,L,15,'like',Sample);
centerIdx = L/2 + 1;

fprintf('Generating A4-style 3D-SIM dataset: zWashout=%.3g um, alphaFloor=%.3g, bgStrength=%.3g\n', ...
    par.zWashout, par.alphaFloor, par.bgStrength);

for z0 = par.zScanRange
    outZ = z0 + centerIdx;
    if outZ < 1 || outZ > L, continue; end
    SampleRel = circshift(temp,z0,3);

    % Optional phase-independent smooth haze object fixed to the sample.
    if par.bgStrength > 0
        BgRel = localSmooth3D(SampleRel, par.bgSigma);
        maxBg = max(BgRel(:));
        if maxBg > 0, BgRel = BgRel ./ maxBg .* max(SampleRel(:)); end
        fBg = fftshift(fftn(ifftshift(BgRel)));
        ImgBg3D = fftshift(ifftn(ifftshift(fBg .* fPSF)));
        bgPlane = real(ImgBg3D(:,:,centerIdx));
    else
        bgPlane = 0;
    end

    if mod(z0,8)==0, fprintf('  z0 = %d\n', z0); end
    for ch = 1:15
        Object = SampleRel .* PatternEff(:,:,:,ch);
        fObject = fftshift(fftn(ifftshift(Object)));
        Img3D = fftshift(ifftn(ifftshift(fObject .* fPSF)));
        rawPlane = real(Img3D(:,:,centerIdx)) + par.bgStrength .* bgPlane;
        ImgSIM(:,:,outZ,ch) = cast(rawPlane,'like',Sample);
    end
end

% Normalize like A4 but keep as single by default for reconstruction.
ImgSIM = ImgSIM ./ max(ImgSIM(:)) .* 65535;
if ~par.useSingle
    ImgSIM = double(ImgSIM);
end

info = struct();
info.PixelSize = par.PixelSize;
info.ccdPxlSize = par.PixelSize;
info.IncidentAngle = par.IncidentAngle;
info.emWavelength = par.emWavelength;
info.zWashout = par.zWashout;
info.alphaFloor = par.alphaFloor;
info.alphaZ = alphaZ;
info.zScanRange = par.zScanRange;
info.bgStrength = par.bgStrength;
info.bgSigma = par.bgSigma;
info.PatternIdealOneAngle = PatternIdealOneAngle;
info.PatternEffOneAngle = PatternEff(:,:,:,1:5);
info.temp = temp;
info.note = 'A4-style 3-beam 3D-SIM dataset with depth-dependent modulation washout.';

if ~isempty(par.saveFile)
    save(par.saveFile,'ImgSIM','info','-v7.3');
    fprintf('Saved %s\n', par.saveFile);
end

if par.showFigures
    figure;
    subplot(221); imagesc(squeeze(ImgSIM(:,:,centerIdx,1))); axis image off; colormap hot; title('Raw XY, ch 1');
    subplot(222); imagesc(squeeze(ImgSIM(:,centerIdx,:,1))'); axis image off; colormap hot; title('Raw XZ, ch 1');
    subplot(223); imagesc(squeeze(alphaZ(:,centerIdx,:))'); axis image off; colormap hot; title('alpha(z) XZ');
    subplot(224); imagesc(squeeze(temp(:,:,centerIdx))); axis image off; colormap hot; title('Sample central XY');
end
end

function B = localSmooth3D(A,sigma)
% Uses imgaussfilt3 if available; otherwise falls back to separable convn.
if exist('imgaussfilt3','file') == 2
    B = imgaussfilt3(A,sigma);
else
    rad = max(1,ceil(3*sigma));
    x = -rad:rad;
    g = exp(-(x.^2)/(2*sigma^2));
    g = g./sum(g);
    B = convn(A,reshape(g,[],1,1),'same');
    B = convn(B,reshape(g,1,[],1),'same');
    B = convn(B,reshape(g,1,1,[]),'same');
end
end
