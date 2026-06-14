% Ye Ma, 2026.06.13
% Department of Biomedical Engineering, Johns Hopkins University


diag = A8_Visualization('Img3DSIM_StrongDegrade.mat');


function diag = A8_Visualization(matFile, opts)
% Visualize, from an A6-generated dataset if available:
%   1) degraded 3D-SIM pattern profile alpha(z)
%   2) 3D sample structure / axial sample distribution
%   3) fluorescence distribution after pattern illumination before PSF/scanning
%   4) final 3D raw image set after scanning/convolution
%   5) IDC1, A1, A2, IDC2, IDC2/IDC1 from the raw phase stack
%
% Optional:
%   opts.angleIndex = 1;
%   opts.phaseIndex = 1;
%   opts.zIndex = [];
%   opts.xIndex = [];
%   opts.yIndex = [];
%   opts.ccdPxlSize = 0.08;
%   opts.recomputePatternIfMissing = true;
%   opts.recomputeIlluminatedVolume = true;
%
% Notes:
%   - This script tries to read variables saved by A6: ImgSIM, info/infoG.
%   - If PatternEff / temp / SampleRel are not saved, it reconstructs the
%     A4-style degraded pattern and sample slab from the parameters.
%   - It does NOT regenerate the whole scanned dataset; it only reconstructs
%     one representative objective-centered illuminated 3D volume for display.

if nargin < 1 || isempty(matFile)
    matFile = 'Img3DSIM_StrongDegrade.mat';
end
if nargin < 2 || isempty(opts)
    opts = struct();
end

opts = setDefault(opts,'angleIndex',1);
opts = setDefault(opts,'phaseIndex',1);
opts = setDefault(opts,'zIndex',[]);
opts = setDefault(opts,'xIndex',[]);
opts = setDefault(opts,'yIndex',[]);
opts = setDefault(opts,'ccdPxlSize',0.08);
opts = setDefault(opts,'IncidentAngle',pi/3);
opts = setDefault(opts,'emWavelength',0.488);
opts = setDefault(opts,'zWashout',0.5);
opts = setDefault(opts,'alphaFloor',0.0);
opts = setDefault(opts,'recomputePatternIfMissing',true);
opts = setDefault(opts,'recomputeIlluminatedVolume',true);
opts = setDefault(opts,'saveDiagMat',true);

S = load(matFile);
if ~isfield(S,'ImgSIM')
    error('MAT file does not contain ImgSIM.');
end
ImgSIM = single(S.ImgSIM);

% Find info struct if present.
info = struct();
fn = fieldnames(S);
for i = 1:numel(fn)
    if isstruct(S.(fn{i}))
        info = S.(fn{i});
        break;
    end
end

% Dimensions: support both [X Y Z 15] and [X Y Z angle phase]
sz = size(ImgSIM);
if numel(sz) == 5
    Lx = sz(1); Ly = sz(2); Lz = sz(3); nAngles = sz(4); nPhases = sz(5);
    getStack = @(ang) squeeze(ImgSIM(:,:,:,ang,:));
elseif numel(sz) == 4
    Lx = sz(1); Ly = sz(2); Lz = sz(3); nChannels = sz(4);
    if nChannels == 15
        nAngles = 3; nPhases = 5;
        getStack = @(ang) ImgSIM(:,:,:,(ang-1)*5 + (1:5));
    else
        nAngles = 1; nPhases = nChannels;
        getStack = @(ang) ImgSIM;
    end
else
    error('Unexpected ImgSIM size. Expected [X Y Z 15] or [X Y Z angle phase].');
end

L = Lx;
if isempty(opts.zIndex), opts.zIndex = round(Lz/2); end
if isempty(opts.xIndex), opts.xIndex = round(Lx/2); end
if isempty(opts.yIndex), opts.yIndex = round(Ly/2); end

ang = min(max(opts.angleIndex,1),nAngles);
ph  = min(max(opts.phaseIndex,1),nPhases);
z0  = opts.zIndex;
x0  = opts.xIndex;
y0  = opts.yIndex;

ccdPxlSize = getInfoOrDefault(info, {'ccdPxlSize','pixelSize','PixelSize'}, opts.ccdPxlSize);
IncidentAngle = getInfoOrDefault(info, {'IncidentAngle','incidentAngle'}, opts.IncidentAngle);
emWavelength = getInfoOrDefault(info, {'emWavelength','wavelength','lambda'}, opts.emWavelength);
zWashout = getInfoOrDefault(info, {'zWashout'}, opts.zWashout);
alphaFloor = getInfoOrDefault(info, {'alphaFloor'}, opts.alphaFloor);

coordX = ((1:Lx)-ceil(Lx/2))*ccdPxlSize;
coordY = ((1:Ly)-ceil(Ly/2))*ccdPxlSize;
z_um   = ((1:Lz)-ceil(Lz/2))*ccdPxlSize;

%% Get or reconstruct sample slab temp
[temp, sampleSource] = getSampleVolume(S, info, Lx, Ly, Lz);

%% Get or reconstruct alphaZ and degraded pattern for one angle / five phases
[PatternIdeal, PatternEff, alphaZ, patternSource] = getPatternVolumes(S, info, Lx, Ly, Lz, nPhases, ang, ccdPxlSize, IncidentAngle, emWavelength, zWashout, alphaFloor, opts.recomputePatternIfMissing);

%% Representative objective-centered sample and illuminated fluorescence volume
% Use scan plane at center z; this matches the objective-centered coordinate.
scanZIndex = z0;
zShift = round(Lz/2) - scanZIndex;
SampleRel = circshift(temp,[0 0 zShift]);
IllumFluo = [];
if opts.recomputeIlluminatedVolume
    IllumFluo = SampleRel .* PatternEff(:,:,:,ph);
end

%% Actual raw scanned image stack for chosen angle
Stack = getStack(ang); % [X Y Z phase]
RawOnePhase = Stack(:,:,:,ph);

%% Lock-in diagnostic terms from actual raw stack
[IDC1, A1, A2, IDC2, ModMap] = computeLockinDiagnostics(Stack);

%% Axial profiles
sampleZ = squeeze(sum(sum(temp,1),2));
sampleZ = normalizeSafe(sampleZ);

alphaZmean = squeeze(mean(mean(alphaZ,1),2));
alphaZmean = normalizeSafe(alphaZmean);

PatternDC = mean(PatternEff,4);
[PAC1, PAC2] = harmonicAmplitudes(PatternEff);
patMod1Z = squeeze(mean(mean(PAC1./(PatternDC+eps('single')),1),2));
patMod2Z = squeeze(mean(mean(PAC2./(PatternDC+eps('single')),1),2));
patMod1Z = normalizeSafe(patMod1Z);
patMod2Z = normalizeSafe(patMod2Z);

if ~isempty(IllumFluo)
    illumZ = squeeze(sum(sum(IllumFluo,1),2));
    illumZ = normalizeSafe(illumZ);
else
    illumZ = nan(size(sampleZ));
end

rawZ = squeeze(mean(mean(RawOnePhase,1),2));
rawZ = normalizeSafe(rawZ);
IDC1z = normalizeSafe(squeeze(mean(mean(IDC1,1),2)));
IDC2z = normalizeSafe(squeeze(mean(mean(IDC2,1),2)));
ModZ  = normalizeSafe(squeeze(mean(mean(ModMap,1),2)));

%% Figure 1: axial comparison
figure('Name','Axial profiles: sample, degraded pattern, illuminated fluorescence, raw scan');
plot(z_um, sampleZ, 'LineWidth', 2); hold on;
plot(z_um, alphaZmean, '--', 'LineWidth', 2);
plot(z_um, patMod1Z, 'LineWidth', 2);
plot(z_um, illumZ, 'LineWidth', 2);
plot(z_um, rawZ, 'LineWidth', 2);
grid on;
xlabel('z (um)'); ylabel('normalized value');
legend('Sample density','alpha(z)','Pattern 1st harmonic contrast','After pattern illumination','Raw scanned image','Location','best');
title(sprintf('Axial comparison, angle %d phase %d',ang,ph));

%% Figure 2: degraded pattern slices
figure('Name','Degraded 3D-SIM pattern');
show3(PatternEff(:,:,:,ph), coordX, coordY, z_um, x0, y0, z0, 'Degraded pattern P_{eff}');

%% Figure 3: sample structure
figure('Name','3D sample structure');
show3(temp, coordX, coordY, z_um, x0, y0, z0, 'Sample volume');

%% Figure 4: fluorescence after pattern illumination before PSF/scanning
if ~isempty(IllumFluo)
    figure('Name','Fluorescence after degraded pattern illumination');
    show3(IllumFluo, coordX, coordY, z_um, x0, y0, z0, 'SampleRel .* P_{eff}');
end

%% Figure 5: actual raw image set after scanning/convolution
figure('Name','Actual raw scanned dataset, one phase');
show3(RawOnePhase, coordX, coordY, z_um, x0, y0, z0, 'Raw scanned image volume');

%% Figure 6: raw phase images at center z
figure('Name','Raw phase images at selected z');
for p = 1:nPhases
    subplot(2,ceil(nPhases/2),p);
    imagesc(coordX,coordY,Stack(:,:,z0,p)); axis image off; colormap hot;
    title(sprintf('Phase %d',p));
end
sgtitle(sprintf('Raw phases at z index %d, angle %d',z0,ang));

%% Figure 7: lock-in diagnostic maps
figure('Name','Lock-in diagnostic maps from actual raw dataset');
subplot(2,3,1); imagesc(coordX,coordY,IDC1(:,:,z0)); axis image off; colormap hot; colorbar; title('IDC1 XY');
subplot(2,3,2); imagesc(coordX,z_um,squeeze(IDC1(:,y0,:))'); axis image; colormap hot; colorbar; title('IDC1 XZ'); xlabel('x (um)'); ylabel('z (um)');
subplot(2,3,3); imagesc(coordY,z_um,squeeze(IDC1(x0,:,:))'); axis image; colormap hot; colorbar; title('IDC1 YZ'); xlabel('y (um)'); ylabel('z (um)');
subplot(2,3,4); imagesc(coordX,coordY,IDC2(:,:,z0)); axis image off; colormap hot; colorbar; title('IDC2 XY');
subplot(2,3,5); imagesc(coordX,z_um,squeeze(IDC2(:,y0,:))'); axis image; colormap hot; colorbar; title('IDC2 XZ'); xlabel('x (um)'); ylabel('z (um)');
subplot(2,3,6); imagesc(coordY,z_um,squeeze(IDC2(x0,:,:))'); axis image; colormap hot; colorbar; title('IDC2 YZ'); xlabel('y (um)'); ylabel('z (um)');

%% Figure 8: ModMap
figure('Name','ModMap = IDC2 / IDC1');
show3(ModMap, coordX, coordY, z_um, x0, y0, z0, 'IDC2 / IDC1');
colorbar;

%% Figure 9: raw-derived axial diagnostics
figure('Name','Raw-derived lock-in axial profiles');
plot(z_um, IDC1z, 'LineWidth', 2); hold on;
plot(z_um, IDC2z, 'LineWidth', 2);
plot(z_um, ModZ, 'LineWidth', 2);
grid on;
xlabel('z (um)'); ylabel('normalized value');
legend('IDC1','IDC2','IDC2/IDC1','Location','best');
title('Lock-in quantities from actual raw phase stack');

%% Output

diag = struct();
diag.matFile = matFile;
diag.info = info;
diag.sampleSource = sampleSource;
diag.patternSource = patternSource;
diag.z_um = z_um;
diag.sampleZ = sampleZ;
diag.alphaZmean = alphaZmean;
diag.patMod1Z = patMod1Z;
diag.patMod2Z = patMod2Z;
diag.illumZ = illumZ;
diag.rawZ = rawZ;
diag.IDC1 = IDC1;
diag.A1 = A1;
diag.A2 = A2;
diag.IDC2 = IDC2;
diag.ModMap = ModMap;
diag.IDC1z = IDC1z;
diag.IDC2z = IDC2z;
diag.ModZ = ModZ;

if opts.saveDiagMat
    save('A8_DatasetVisualization_Diagnostics.mat','diag','-v7.3');
    fprintf('Saved A8_DatasetVisualization_Diagnostics.mat\n');
end

fprintf('Sample source: %s\n', sampleSource);
fprintf('Pattern source: %s\n', patternSource);
end

%% Helper functions
function opts = setDefault(opts,name,value)
if ~isfield(opts,name) || isempty(opts.(name))
    opts.(name) = value;
end
end

function val = getInfoOrDefault(info, names, defaultVal)
val = defaultVal;
for i = 1:numel(names)
    if isfield(info,names{i})
        val = info.(names{i});
        return;
    end
end
end

function [temp, source] = getSampleVolume(S, info, Lx, Ly, Lz)
source = 'reconstructed from Sample1000Tubule_3D_256Pixels.mat or zeros';

candidateNames = {'temp','SampleTemp','SampleSlab','sampleVolume','SampleVolume','SampleRel0','Sample'};
for i = 1:numel(candidateNames)
    nm = candidateNames{i};
    if isfield(S,nm)
        A = single(S.(nm));
        if isequal(size(A),[Lx Ly Lz])
            temp = A; source = ['loaded S.' nm]; return;
        end
    end
    if isfield(info,nm)
        A = single(info.(nm));
        if isequal(size(A),[Lx Ly Lz])
            temp = A; source = ['loaded info.' nm]; return;
        end
    end
end

if exist('Sample1000Tubule_3D_256Pixels.mat','file')
    T = load('Sample1000Tubule_3D_256Pixels.mat');
    if isfield(T,'Sample')
        Sample = single(T.Sample);
        temp = zeros(Lx,Ly,Lz,'single');
        z1 = max(1,round(Lz/4)+1);
        z2 = min(Lz,round(3*Lz/4));
        if isequal(size(Sample),[Lx Ly Lz])
            temp(:,:,z1:z2) = Sample(:,:,z1:z2);
        else
            sx = min(Lx,size(Sample,1)); sy = min(Ly,size(Sample,2)); sz = min(Lz,size(Sample,3));
            temp(1:sx,1:sy,1:sz) = Sample(1:sx,1:sy,1:sz);
        end
        source = 'loaded Sample1000Tubule_3D_256Pixels.mat and embedded central slab';
        return;
    end
end

temp = zeros(Lx,Ly,Lz,'single');
warning('Could not find sample volume. Using zeros.');
end

function [PatternIdeal, PatternEff, alphaZ, source] = getPatternVolumes(S, info, Lx, Ly, Lz, nPhases, angleIndex, ccdPxlSize, IncidentAngle, emWavelength, zWashout, alphaFloor, recompute)
source = 'recomputed A4-style degraded pattern';
PatternIdeal = [];
PatternEff = [];
alphaZ = [];

% Try saved fields first.
for nm = {'PatternEff','Pattern3DSIM_Eff','PatternDegraded'}
    name = nm{1};
    if isfield(S,name), A = single(S.(name)); else, A = []; end
    if isempty(A) && isfield(info,name), A = single(info.(name)); end
    if ~isempty(A)
        if ndims(A)==4 && size(A,4)>=nPhases
            PatternEff = A(:,:,:,1:nPhases); source = ['loaded ' name];
        elseif ndims(A)==5 && size(A,4)>=angleIndex && size(A,5)>=nPhases
            PatternEff = squeeze(A(:,:,:,angleIndex,1:nPhases)); source = ['loaded ' name];
        end
    end
end

for nm = {'PatternIdeal','Pattern3DSIM','Pattern'}
    name = nm{1};
    if isfield(S,name), A = single(S.(name)); else, A = []; end
    if isempty(A) && isfield(info,name), A = single(info.(name)); end
    if ~isempty(A)
        if ndims(A)==4 && size(A,4)>=nPhases
            PatternIdeal = A(:,:,:,1:nPhases);
        elseif ndims(A)==5 && size(A,4)>=angleIndex && size(A,5)>=nPhases
            PatternIdeal = squeeze(A(:,:,:,angleIndex,1:nPhases));
        end
    end
end

for nm = {'alphaZ','alpha','modulationAlpha'}
    name = nm{1};
    if isfield(S,name), A = single(S.(name)); else, A = []; end
    if isempty(A) && isfield(info,name), A = single(info.(name)); end
    if ~isempty(A)
        alphaZ = A;
        if isvector(alphaZ)
            alphaZ = reshape(alphaZ,1,1,[]);
            alphaZ = repmat(alphaZ,[Lx Ly 1]);
        end
    end
end

if ~isempty(PatternEff) && isempty(PatternIdeal)
    PatternIdeal = PatternEff;
end

if ~isempty(PatternEff) && isempty(alphaZ)
    alphaZ = ones(Lx,Ly,Lz,'single');
end

if ~isempty(PatternEff) && ~recompute
    return;
end

if isempty(PatternEff)
    coord = ((1:Lx)-ceil(Lx/2))*ccdPxlSize;
    [xxx,yyy,zzz] = meshgrid(coord,coord,coord);
    k = 2*pi/emWavelength;
    kr = k*sin(IncidentAngle);
    kz = k*cos(IncidentAngle);
    thetaSet = [0, pi/3, 2*pi/3];
    theta = thetaSet(min(angleIndex,numel(thetaSet)));
    kx = kr*cos(theta); ky = kr*sin(theta);
    phiVec = 2*pi/nPhases*(0:nPhases-1);
    alphaZ = alphaFloor + (1-alphaFloor).*exp(-(zzz./zWashout).^2);
    alphaZ = single(alphaZ);
    PatternIdeal = zeros(Lx,Ly,Lz,nPhases,'single');
    PatternEff = zeros(Lx,Ly,Lz,nPhases,'single');
    for p = 1:nPhases
        phi1 = phiVec(p);
        phi2 = phi1/2;
        E1 = exp(1j.*(kx.*xxx + ky.*yyy + phi1) + 1j.*kz.*zzz);
        E2 = exp(1j.*(k.*zzz + phi2));
        E3 = exp(-1j.*(kx.*xxx + ky.*yyy) + 1j.*kz.*zzz);
        P = abs(E1+E2+E3).^2;
        P = P ./ mean(P(:));
        PatternIdeal(:,:,:,p) = single(P);
        PatternEff(:,:,:,p) = single(1 + alphaZ.*(P-1));
    end
end
end

function [IDC1, A1, A2, IDC2, ModMap] = computeLockinDiagnostics(Stack)
Stack = single(Stack);
nPhases = size(Stack,4);
phi = 2*pi*(0:nPhases-1)/nPhases;
IDC1 = mean(Stack,4);
C1 = zeros(size(IDC1),'single'); S1 = zeros(size(IDC1),'single');
C2 = zeros(size(IDC1),'single'); S2 = zeros(size(IDC1),'single');
for p = 1:nPhases
    I = Stack(:,:,:,p);
    C1 = C1 + I*cos(phi(p));
    S1 = S1 + I*sin(phi(p));
    C2 = C2 + I*cos(2*phi(p));
    S2 = S2 + I*sin(2*phi(p));
end
A1 = (2/nPhases)*sqrt(C1.^2 + S1.^2);
A2 = (2/nPhases)*sqrt(C2.^2 + S2.^2);
IDC2 = A1 + A2;
ModMap = IDC2 ./ (IDC1 + eps('single'));
end

function [A1,A2] = harmonicAmplitudes(P)
nPhases = size(P,4);
phi = 2*pi*(0:nPhases-1)/nPhases;
DC = mean(P,4);
C1 = zeros(size(DC),'single'); S1 = zeros(size(DC),'single');
C2 = zeros(size(DC),'single'); S2 = zeros(size(DC),'single');
for p = 1:nPhases
    I = single(P(:,:,:,p));
    C1 = C1 + I*cos(phi(p));
    S1 = S1 + I*sin(phi(p));
    C2 = C2 + I*cos(2*phi(p));
    S2 = S2 + I*sin(2*phi(p));
end
A1 = (2/nPhases)*sqrt(C1.^2 + S1.^2);
A2 = (2/nPhases)*sqrt(C2.^2 + S2.^2);
end

function y = normalizeSafe(x)
x = double(squeeze(x));
mx = max(x(:));
if mx > 0
    y = x ./ mx;
else
    y = x;
end
end

function show3(V, coordX, coordY, z_um, x0, y0, z0, ttl)
V = squeeze(V);
subplot(1,3,1);
imagesc(coordX,coordY,V(:,:,z0)); axis image off; colormap hot; colorbar;
title([ttl ' XY']);
subplot(1,3,2);
imagesc(coordX,z_um,squeeze(V(:,y0,:))'); axis image; colormap hot; colorbar;
xlabel('x (um)'); ylabel('z (um)'); title([ttl ' XZ']);
subplot(1,3,3);
imagesc(coordY,z_um,squeeze(V(x0,:,:))'); axis image; colormap hot; colorbar;
xlabel('y (um)'); ylabel('z (um)'); title([ttl ' YZ']);
end
