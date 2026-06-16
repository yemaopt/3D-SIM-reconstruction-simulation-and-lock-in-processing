% Ye Ma, 2026.06.16
% Compare standard 3D-SIM, notch 3D-SIM, and two Lock-in 3D-SIM variants:
%   1) harmonic-extraction IDC2
%   2) demix-order IDC2 scaled to harmonic-equivalent amplitude
%
% Outputs:
%   Recon.standard_no_notch
%   Recon.standard_notch
%   Recon.lockin_harmonic_no_notch
%   Recon.lockin_demix_no_notch
%
% Notes:
%   - harmonic IDC2 uses direct phase harmonic extraction with the 2/5 factor.
%   - demix IDC2 uses A5-style demixing, then multiplies by sqrt(2) so that
%     it matches the harmonic amplitude scale in the ideal case.

A7_Compare_Standard3DSIMRecon_NoNotch_Notch_Lockin3DSIM_func();

function [Recon, info] = A7_Compare_Standard3DSIMRecon_NoNotch_Notch_Lockin3DSIM_func(imgFile, varargin)

if nargin < 1 || isempty(imgFile)
    imgFile = 'Img3DSIM_StrongDegrade.mat';
end

p = inputParser;
p.addParameter('epsilonLock',1,@isnumeric);
p.addParameter('gammaLock',1.0,@isnumeric);
p.addParameter('lockinUseSecondOrder',true,@islogical);
p.addParameter('clipNegative',true,@islogical);
p.addParameter('notchSigma',15,@isnumeric);
p.addParameter('notchA',0.9,@isnumeric);
p.addParameter('PixelSizeAfterInterp',0.08/2,@isnumeric);
p.addParameter('IncidentAngle',pi/3,@isnumeric);
p.addParameter('emWavelength',0.488,@isnumeric);
p.addParameter('showFigures',true,@islogical);
p.addParameter('saveFile','A7_compare_lockin_harmonic_demix_v2_results.mat',@ischar);
p.parse(varargin{:});
par = p.Results;

S = load(imgFile);
if ~isfield(S,'ImgSIM'), error('Input file must contain ImgSIM'); end
ImgSIM0 = single(S.ImgSIM);

% Use dataset metadata if available.
if isfield(S,'info')
    if isfield(S.info,'IncidentAngle'), par.IncidentAngle = S.info.IncidentAngle; end
    if isfield(S.info,'emWavelength'), par.emWavelength = S.info.emWavelength; end
    if isfield(S.info,'PixelSize'), par.PixelSizeAfterInterp = S.info.PixelSize/2; end
end

fprintf('\nRunning reconstructions on %s\n', imgFile);
Recon = struct();

Recon.standard_no_notch = reconstructA5Style(ImgSIM0, false, par, 'standard_no_notch');
Recon.standard_notch    = reconstructA5Style(ImgSIM0, true,  par, 'standard_notch');

[ImgSIM_lock_harm, diagH] = apply3DLockinRaw_harmonic(ImgSIM0, par);
Recon.lockin_harmonic_no_notch = reconstructA5Style(ImgSIM_lock_harm, false, par, 'lockin_harmonic_no_notch');

[ImgSIM_lock_demix, diagD] = apply3DLockinRaw_demix(ImgSIM0, par);
Recon.lockin_demix_no_notch = reconstructA5Style(ImgSIM_lock_demix, false, par, 'lockin_demix_no_notch');

info = struct();
info.imgFile = imgFile;
info.epsilonLock = par.epsilonLock;
info.gammaLock = par.gammaLock;
info.lockinUseSecondOrder = par.lockinUseSecondOrder;
info.notchSigma = par.notchSigma;
info.notchA = par.notchA;
info.lockin_harmonic_diag = diagH;
info.lockin_demix_diag = diagD;
info.note = ['Comparison of standard no-notch, standard notch, ', ...
             'raw Lock-in harmonic IDC2 no-notch, and raw Lock-in demix IDC2 no-notch.'];
if isfield(S,'info'), info.datasetInfo = S.info; end

if ~isempty(par.saveFile)
    save(par.saveFile,'Recon','info','-v7.3');
    fprintf('Saved %s\n', par.saveFile);
end

if par.showFigures
    showComparison(Recon);
    showLockinDiagnostics(diagH, diagD);
end
end

function ReconTrunc = reconstructA5Style(ImgSIM, useNotch, par, label)
thetaSet = [0, 1/3, 2/3].*pi;
phiVec = 2*pi/5.*(0:4);
PtnNum = 5;

[fftDim_LR_XY,~,fftDim_Z,~] = size(ImgSIM);
fftDim_HR_XY = fftDim_LR_XY * 2;
LRinHRRange_XY = (-fftDim_LR_XY/2:fftDim_LR_XY/2-1) + (fftDim_HR_XY/2+1);

[XXX,YYY,~] = meshgrid(-fftDim_HR_XY/2:fftDim_HR_XY/2-1, ...
                       -fftDim_HR_XY/2:fftDim_HR_XY/2-1, ...
                       -fftDim_Z/2:fftDim_Z/2-1);
XXX = single(XXX); YYY = single(YYY);

if useNotch
    NotchFilter = single(1 - par.notchA .* exp(-(XXX.^2+YYY.^2)./(2*par.notchSigma^2)));
else
    NotchFilter = single(ones(fftDim_HR_XY,fftDim_HR_XY,fftDim_Z));
end

k = 2*pi/par.emWavelength;
kr = k*sin(par.IncidentAngle);

ModDepth = 1;
M = diag([1,1/3*ModDepth,1/3*ModDepth,1/3*ModDepth,1/3*ModDepth]) * ...
    exp([0; -0.5j; 0.5j; -1j; 1j] * phiVec);
[U,Sv,V] = svd(M);
Minv = single(V * (Sv\eye(5)) * U');

ReconSum = complex(zeros(fftDim_HR_XY,fftDim_HR_XY,fftDim_Z,'single'));
fprintf('  reconstructing %s, useNotch=%d\n', label, useNotch);

for thetaIndx = 1:3
    theta = thetaSet(thetaIndx);
    kx = kr*cos(theta);
    ky = kr*sin(theta);

    phaseF = complex(zeros(fftDim_LR_XY^2*fftDim_Z,PtnNum,'single'));
    for p = 1:PtnNum
        tmp = ImgSIM(:,:,:,(thetaIndx-1)*5+p);
        F = fftshift(fftn(ifftshift(tmp)));
        phaseF(:,p) = F(:);
    end

    sol = phaseF * Minv;
    bandsLR = reshape(sol,fftDim_LR_XY,fftDim_LR_XY,fftDim_Z,PtnNum);
    clear phaseF sol

    for b = 1:5
        if thetaIndx > 1 && b == 1
            continue;
        end

        tmpFT = complex(zeros(fftDim_HR_XY,fftDim_HR_XY,fftDim_Z,'single'));
        tmpFT(LRinHRRange_XY,LRinHRRange_XY,:) = bandsLR(:,:,:,b);

        % A5-style notch: zero and +/-2 bands.
        if useNotch && (b == 1 || b == 4 || b == 5)
            tmpFT = tmpFT .* NotchFilter;
        end

        tmpImg = fftshift(ifftn(ifftshift(tmpFT)));
        switch b
            case 1
                ramp = 1;
            case 2
                ramp = exp(1j*(+kx.*XXX.*par.PixelSizeAfterInterp + ky.*YYY.*par.PixelSizeAfterInterp));
            case 3
                ramp = exp(1j*(-kx.*XXX.*par.PixelSizeAfterInterp - ky.*YYY.*par.PixelSizeAfterInterp));
            case 4
                ramp = exp(1j*(+2*kx.*XXX.*par.PixelSizeAfterInterp + 2*ky.*YYY.*par.PixelSizeAfterInterp));
            case 5
                ramp = exp(1j*(-2*kx.*XXX.*par.PixelSizeAfterInterp - 2*ky.*YYY.*par.PixelSizeAfterInterp));
        end
        ReconSum = ReconSum + tmpImg .* ramp;
    end
end

ReconAbs = abs(ReconSum);
z1 = fftDim_Z/4 + 1;
z2 = 3*fftDim_Z/4;
ReconTrunc = ReconAbs(:,:,z1:z2);
end

function [ImgSIM_lock, diag1] = apply3DLockinRaw_harmonic(ImgSIM, par)
% Direct lock-in harmonic extraction from raw phase stack.
% This is closest to the equations:
%   IDC1 = mean(Ij)
%   A1 = (2/5)*sqrt((sum Ij cos(phi/2))^2 + (sum Ij sin(phi/2))^2)
%   A2 = (2/5)*sqrt((sum Ij cos(phi))^2   + (sum Ij sin(phi))^2)
%   IDC2 = A1 + A2
phiVec = 2*pi/5.*(0:4);
PtnNum = 5;
[Ny,Nx,Nz,Nc] = size(ImgSIM);
if Nc ~= 15, error('ImgSIM must have 15 channels = 3 angles x 5 phases'); end

ImgSIM_lock = zeros(size(ImgSIM),'like',ImgSIM);
diag1 = struct();
fprintf('  applying raw 3D Lock-in HARMONIC: epsilon=%.3g, gamma=%.3g\n', par.epsilonLock, par.gammaLock);

for thetaIndx = 1:3
    raw = ImgSIM(:,:,:,(thetaIndx-1)*5+(1:5));
    IDC1 = mean(raw,4);

    C1 = zeros(Ny,Nx,Nz,'like',ImgSIM);
    S1 = zeros(Ny,Nx,Nz,'like',ImgSIM);
    C2 = zeros(Ny,Nx,Nz,'like',ImgSIM);
    S2 = zeros(Ny,Nx,Nz,'like',ImgSIM);

    for p = 1:PtnNum
        ph = phiVec(p).*2;
        I = raw(:,:,:,p);
        C1 = C1 + I .* cos(ph/2);
        S1 = S1 + I .* sin(ph/2);
        C2 = C2 + I .* cos(ph);
        S2 = S2 + I .* sin(ph);
    end

    A1 = (2/PtnNum) .* sqrt(C1.^2 + S1.^2);
    if par.lockinUseSecondOrder
        A2 = (2/PtnNum) .* sqrt(C2.^2 + S2.^2);
    else
        A2 = zeros(size(A1),'like',A1);
    end
    IDC2 = single(A1 + A2);

    for p = 1:PtnNum
        out = raw(:,:,:,p) - par.epsilonLock .* IDC1 + par.gammaLock .* IDC2;
        if par.clipNegative, out(out < 0) = 0; end
        ImgSIM_lock(:,:,:,(thetaIndx-1)*5+p) = out;
    end

    if thetaIndx == 1
        diag1.IDC1 = IDC1;
        diag1.A1 = single(A1);
        diag1.A2 = single(A2);
        diag1.IDC2 = IDC2;
        diag1.modulationRatio = IDC2 ./ (IDC1 + eps('single'));
    end
end
end

function [ImgSIM_lock, diag1] = apply3DLockinRaw_demix(ImgSIM, par)
% Demix-order IDC2. In the ideal case, harmonic amplitude = sqrt(2)*demix amplitude.
% Therefore IDC2_demix_scaled = sqrt(2)*(sqrt(|B+1|^2+|B-1|^2)+sqrt(|B+2|^2+|B-2|^2)).
% Do NOT multiply by 2/5 here because Minv already demixes the phase stack.
phiVec = 2*pi/5.*(0:4);
PtnNum = 5;
[Ny,Nx,Nz,Nc] = size(ImgSIM);
if Nc ~= 15, error('ImgSIM must have 15 channels = 3 angles x 5 phases'); end

ModDepth = 1;
M = diag([1,1/3*ModDepth,1/3*ModDepth,1/3*ModDepth,1/3*ModDepth]) * ...
    exp([0; -0.5j; 0.5j; -1j; 1j] * phiVec);
[U,Sv,V] = svd(M);
Minv = single(V * (Sv\eye(5)) * U');

ImgSIM_lock = zeros(size(ImgSIM),'like',ImgSIM);
diag1 = struct();
fprintf('  applying raw 3D Lock-in DEMIX: epsilon=%.3g, gamma=%.3g\n', par.epsilonLock, par.gammaLock);

for thetaIndx = 1:3
    raw = ImgSIM(:,:,:,(thetaIndx-1)*5+(1:5));
    IDC1 = mean(raw,4);

    R = reshape(raw,[],PtnNum);
    B = complex(R) * Minv;
    B = reshape(B,Ny,Nx,Nz,PtnNum);

    A1_demix = sqrt(abs(B(:,:,:,2)).^2 + abs(B(:,:,:,3)).^2);
    if par.lockinUseSecondOrder
        A2_demix = sqrt(abs(B(:,:,:,4)).^2 + abs(B(:,:,:,5)).^2);
    else
        A2_demix = zeros(size(A1_demix),'like',A1_demix);
    end

    A1 = sqrt(2)./3 .* A1_demix;
    A2 = sqrt(2)./3 .* A2_demix;
    IDC2 = single(A1 + A2);

    for p = 1:PtnNum
        out = raw(:,:,:,p) - par.epsilonLock .* IDC1 + par.gammaLock .* IDC2;
        if par.clipNegative, out(out < 0) = 0; end
        ImgSIM_lock(:,:,:,(thetaIndx-1)*5+p) = out;
    end

    if thetaIndx == 1
        diag1.IDC1 = IDC1;
        diag1.A1 = single(A1);
        diag1.A2 = single(A2);
        diag1.IDC2 = IDC2;
        diag1.modulationRatio = IDC2 ./ (IDC1 + eps('single'));
    end
end
end

function showComparison(Recon)
fields = {'standard_no_notch','standard_notch','lockin_harmonic_no_notch','lockin_demix_no_notch'};
titles = {'Standard, no notch','Standard + notch','Lock-in harmonic, no notch','Lock-in demix, no notch'};
figure('Name','Standard vs notch vs two Lock-in IDC2 definitions');
for i = 1:numel(fields)
    R = Recon.(fields{i});
    cy = round(size(R,1)/2); cx = round(size(R,2)/2); cz = round(size(R,3)/2);
    subplot(numel(fields),3,(i-1)*3+1); imagesc(R(:,:,cz)); axis image off; colormap hot; title([titles{i} ' XY']);
    subplot(numel(fields),3,(i-1)*3+2); imagesc(squeeze(R(:,cx,:))'); axis image off; colormap hot; title('XZ');
    subplot(numel(fields),3,(i-1)*3+3); imagesc(squeeze(R(cy,:,:))'); axis image off; colormap hot; title('YZ');
end
end

function showLockinDiagnostics(H,D)
z = round(size(H.IDC1,3)/2);
figure('Name','Lock-in IDC2 diagnostics: harmonic vs demix');
subplot(2,4,1); imagesc(H.IDC1(:,:,z)); axis image off; colormap hot; colorbar; title('IDC1');
subplot(2,4,2); imagesc(H.IDC2(:,:,z)); axis image off; colormap hot; colorbar; title('IDC2 harmonic');
subplot(2,4,3); imagesc(D.IDC2(:,:,z)); axis image off; colormap hot; colorbar; title('IDC2 demix scaled');
subplot(2,4,4); imagesc(abs(H.IDC2(:,:,z)-D.IDC2(:,:,z))); axis image off; colormap hot; colorbar; title('|difference|');
subplot(2,4,5); imagesc(H.modulationRatio(:,:,z)); axis image off; colormap hot; colorbar; title('Harmonic IDC2/IDC1');
subplot(2,4,6); imagesc(D.modulationRatio(:,:,z)); axis image off; colormap hot; colorbar; title('Demix IDC2/IDC1');
subplot(2,4,7:8);
plot(squeeze(mean(mean(H.modulationRatio,1),2)),'LineWidth',2); hold on;
plot(squeeze(mean(mean(D.modulationRatio,1),2)),'LineWidth',2);
grid on; legend('harmonic','demix scaled'); xlabel('z slice'); ylabel('mean IDC2/IDC1'); title('Mean modulation ratio vs z');
end
