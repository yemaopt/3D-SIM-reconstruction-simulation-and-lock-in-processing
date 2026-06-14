% Ye Ma, 2026.06.13
% Department of Biomedical Engineering, Johns Hopkins University

[Recon, info] = A7_Compare_Standard3DSIMRecon_NoNotch_Notch_Lockin3DSIM_func('Img3DSIM_StrongDegrade.mat');

function [Recon, info] = A7_Compare_Standard3DSIMRecon_NoNotch_Notch_Lockin3DSIM_func(imgFile, varargin)
% Compare three A5-style 3D-SIM reconstructions from the same raw dataset:
%   1) standard recovery without notch filter
%   2) standard recovery with A5 notch filter
%   3) 3D Lock-in-SIM raw filtering + standard recovery without notch
%
% The 3D Lock-in raw filter follows the Lock-in-SIM spirit:
%   Ilock_j = I_j - epsilon * IDC1 + gamma * IDC2
% where IDC1 is the phase average, and IDC2 is a modulation-amplitude
% estimate computed from the 5-phase 3D-SIM demodulation model.
%
% Example:
%   [Recon,info] = A7_Compare_Standard3DSIMRecon_NoNotch_Notch_Lockin3DSIM_func('Img3DSIM_StrongDegrade.mat');
%
% Outputs:
%   Recon.standard_no_notch
%   Recon.standard_notch
%   Recon.lockin_no_notch
%

if nargin < 1 || isempty(imgFile), imgFile = 'Img3DSIM_StrongDegrade.mat'; end

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
p.addParameter('saveFile','A6_compare_v5_results.mat',@ischar);
p.parse(varargin{:});
par = p.Results;

S = load(imgFile);
if ~isfield(S,'ImgSIM'), error('Input file must contain ImgSIM'); end
ImgSIM0 = single(S.ImgSIM);

fprintf('\nRunning three reconstructions on %s\n', imgFile);
Recon = struct();

Recon.standard_no_notch = reconstructA5Style(ImgSIM0, false, par, 'standard_no_notch');
Recon.standard_notch    = reconstructA5Style(ImgSIM0, true,  par, 'standard_notch');

ImgSIM_lock = apply3DLockinRaw(ImgSIM0, par);
Recon.lockin_no_notch   = reconstructA5Style(ImgSIM_lock, false, par, 'lockin_no_notch');

info = struct();
info.imgFile = imgFile;
info.epsilonLock = par.epsilonLock;
info.gammaLock = par.gammaLock;
info.lockinUseSecondOrder = par.lockinUseSecondOrder;
info.notchSigma = par.notchSigma;
info.notchA = par.notchA;
info.note = ['Comparison of standard no-notch, standard with A5 notch, ', ...
             'and raw Lock-in 3D-SIM no-notch.'];
if isfield(S,'info'), info.datasetInfo = S.info; end

if ~isempty(par.saveFile)
    save(par.saveFile,'Recon','info','-v7.3');
    fprintf('Saved %s\n', par.saveFile);
end

if par.showFigures
    showComparison(Recon);
end
end

function ReconTrunc = reconstructA5Style(ImgSIM, useNotch, par, label)
% Memory-light A5-style recon: demix, HR insert, optional notch, phase-ramp,
% sum selected bands [1:5,7:10,12:15].
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
    NotchFilter = 1 - par.notchA .* exp(-(XXX.^2+YYY.^2)./(2*par.notchSigma^2));
    NotchFilter = single(NotchFilter);
else
    NotchFilter = single(ones(fftDim_HR_XY,fftDim_HR_XY,fftDim_Z));
end

k = 2*pi/par.emWavelength;
kr = k*sin(par.IncidentAngle);

% Same demodulation matrix as A5.
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
        % Skip repeated zero orders for theta 2 and 3, as in A5.
        if thetaIndx > 1 && b == 1
            continue;
        end

        tmpFT = complex(zeros(fftDim_HR_XY,fftDim_HR_XY,fftDim_Z,'single'));
        tmpFT(LRinHRRange_XY,LRinHRRange_XY,:) = bandsLR(:,:,:,b);

        % A5 applies notch to 0 and +/-2 bands only.
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

function ImgSIM_lock = apply3DLockinRaw(ImgSIM, par)
% 3D Lock-in raw filtering. Uses the same 5-phase demodulation matrix as A5
% in image space to estimate modulation amplitude IDC2.
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
fprintf('  applying raw 3D Lock-in: epsilon=%.3g, gamma=%.3g\n', par.epsilonLock, par.gammaLock);

for thetaIndx = 1:3
    raw = ImgSIM(:,:,:,(thetaIndx-1)*5+(1:5));
    IDC1 = mean(raw,4);

    % Per-voxel phase demodulation into 0,+1,-1,+2,-2 coefficients.
    R = reshape(raw,[],PtnNum);
    B = complex(R) * Minv;
    B = reshape(B,Ny,Nx,Nz,PtnNum);

    A1 = (2/5).*(sqrt(abs(B(:,:,:,2)).^2 + abs(B(:,:,:,3)).^2));
    if par.lockinUseSecondOrder
        A2 = (2/5).*(sqrt(abs(B(:,:,:,4)).^2 + abs(B(:,:,:,5)).^2));
    else
        A2 = 0;
    end
    IDC2 = single(A1 + A2);

    for p = 1:5
        out = raw(:,:,:,p) - par.epsilonLock .* IDC1 + par.gammaLock .* IDC2;
        if par.clipNegative
            out(out < 0) = 0;
        end
        ImgSIM_lock(:,:,:,(thetaIndx-1)*5+p) = out;
    end
end
end

function showComparison(Recon)
fields = {'standard_no_notch','standard_notch','lockin_no_notch'};
titles = {'Standard, no notch','Standard + notch','Lock-in raw, no notch'};
figure;
for i = 1:3
    R = Recon.(fields{i});
    cy = round(size(R,1)/2); cx = round(size(R,2)/2); cz = round(size(R,3)/2);
    % Use independent contrast per method for morphology comparison.
    subplot(3,3,(i-1)*3+1); imagesc(R(:,:,cz)); axis image off; colormap hot; title([titles{i} ' XY']);
    subplot(3,3,(i-1)*3+2); imagesc(squeeze(R(:,cx,:))'); axis image off; colormap hot; title('XZ');
    subplot(3,3,(i-1)*3+3); imagesc(squeeze(R(cy,:,:))'); axis image off; colormap hot; title('YZ');
end
end
