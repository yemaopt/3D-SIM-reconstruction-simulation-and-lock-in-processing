% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab, 2020.08.20
% Department of Biomedical Engineering, Johns Hopkins University


clear
clc
%close all

load Img3DSIM.mat

thetaSet = [0, 1/3, 2/3].*pi;
phiVec = 2*pi/5.*(0:1:4);

[fftDim_LR_XY,fftDim_Z,~]=size(ImgSIM);
fftDim_HR_XY = fftDim_LR_XY * 2;
LRinHRRange_XY = (-fftDim_LR_XY/2:1:fftDim_LR_XY/2-1)+(fftDim_HR_XY/2+1);

PtnNum = 5;
FT_ImgSIM_HR_SeparatedBand = zeros(fftDim_HR_XY, fftDim_HR_XY, fftDim_Z, PtnNum * length(thetaSet));
fftn_tmpImg_Vectorized = zeros( fftDim_LR_XY^2*fftDim_Z, PtnNum );
for thetaIndx = 1:1:3
    thisTheta = thetaSet(thetaIndx);
    for PhaseIndx = 1:1:PtnNum
        tmpImg = ImgSIM(:,:,:, (thetaIndx-1)*5 + PhaseIndx);
        fftn_tmpImg = fftshift(fftn( ifftshift( tmpImg ) ));
        fftn_tmpImg_Vectorized(:, PhaseIndx) = fftn_tmpImg(:);
    end
    ModDepth=1;
    %The order of the solution:
    %0,+1,-1,+2,-2 order
    rmTransMtx = diag([1,1/3*ModDepth,1/3*ModDepth,1/3*ModDepth,1/3*ModDepth])*exp( [0; -0.5j; 0.5j; -1j; 1j] * phiVec );
    [rmU, rmS, rmV] = svd(rmTransMtx); % rmS will be 5x5 square
    rmTransMtx_pseudoInv = rmV * (rmS\eye(5)) * (rmU');
    fftn_Solution = fftn_tmpImg_Vectorized * rmTransMtx_pseudoInv;
    FT_ImgSIM_HR_SeparatedBand(LRinHRRange_XY, LRinHRRange_XY, : ,(thetaIndx-1)*5+(1:5)) = reshape(fftn_Solution, fftDim_LR_XY, fftDim_LR_XY, fftDim_Z, PtnNum);
end

figure;ii=12;
subplot(131);imagesc(log(abs(squeeze(FT_ImgSIM_HR_SeparatedBand(:,:,128,ii)))+1));colormap(hot);axis equal;
subplot(132);imagesc(log(abs(squeeze(FT_ImgSIM_HR_SeparatedBand(:,256,:,ii))')+1));colormap(hot);axis equal;
subplot(133);imagesc(log(abs(squeeze(FT_ImgSIM_HR_SeparatedBand(256,:,:,ii))')+1));colormap(hot);axis equal;

a=log(abs(squeeze(FT_ImgSIM_HR_SeparatedBand(:,256,:,ii))')+1);

IncidentAngle = pi/3;
ccdPxlSize_AfterInterp = 0.08/2; % in um
emWavelength = 0.488; % in um
k = 2*pi/emWavelength;
kr = k*sin(IncidentAngle);

[XXX, YYY, ZZZ] = meshgrid( -fftDim_HR_XY/2:1:fftDim_HR_XY/2-1, -fftDim_HR_XY/2:1:fftDim_HR_XY/2-1, -fftDim_Z/2:1:fftDim_Z/2-1  );
ImgSIM_HR_SeparatedBand = zeros(fftDim_HR_XY, fftDim_HR_XY, fftDim_Z, PtnNum * length(thetaSet));


%% Notch filter to suppress the high intensity at the origin of the -2 +2 order to make the sum spectrum smooth and free of background artifact
sigma=15;A=0.9;
NotchFilter=1-A.*exp(-(XXX.^2+YYY.^2)./2./sigma^2);
%%
for thetaIndx = 1:3
    theta = thetaSet(thetaIndx);
    kx = kr*cos(theta);
    ky = kr*sin(theta);
    ImgSIM_HR_SeparatedBand(:,:,:,1+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( NotchFilter.*FT_ImgSIM_HR_SeparatedBand(:,:,:,1+5*(thetaIndx-1)))) );   
    ImgSIM_HR_SeparatedBand(:,:,:,2+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( FT_ImgSIM_HR_SeparatedBand(:,:,:,2+5*(thetaIndx-1)))) ).* exp(1j * (+kx.*XXX.*ccdPxlSize_AfterInterp + ky.*YYY.*ccdPxlSize_AfterInterp));
    ImgSIM_HR_SeparatedBand(:,:,:,3+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( FT_ImgSIM_HR_SeparatedBand(:,:,:,3+5*(thetaIndx-1)))) ).* exp(1j * (-kx.*XXX.*ccdPxlSize_AfterInterp - ky.*YYY.*ccdPxlSize_AfterInterp));
    ImgSIM_HR_SeparatedBand(:,:,:,4+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( NotchFilter.*FT_ImgSIM_HR_SeparatedBand(:,:,:,4+5*(thetaIndx-1)))) ).* exp(1j * (+2.*kx.*XXX.*ccdPxlSize_AfterInterp + 2.*ky.*YYY.*ccdPxlSize_AfterInterp));
    ImgSIM_HR_SeparatedBand(:,:,:,5+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( NotchFilter.*FT_ImgSIM_HR_SeparatedBand(:,:,:,5+5*(thetaIndx-1)))) ).* exp(1j * (-2.*kx.*XXX.*ccdPxlSize_AfterInterp - 2.*ky.*YYY.*ccdPxlSize_AfterInterp));
%     ImgSIM_HR_SeparatedBand(:,:,:,4+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( FT_ImgSIM_HR_SeparateBand(:,:,:,4+5*(thetaIndx-1)))) ).* exp(1j * (+2.*kx.*XXX.*ccdPxlSize_AfterInterp + 2.*ky.*YYY.*ccdPxlSize_AfterInterp));
%     ImgSIM_HR_SeparatedBand(:,:,:,5+5*(thetaIndx-1))=fftshift( ifftn( ifftshift( FT_ImgSIM_HR_SeparateBand(:,:,:,5+5*(thetaIndx-1)))) ).* exp(1j * (-2.*kx.*XXX.*ccdPxlSize_AfterInterp - 2.*ky.*YYY.*ccdPxlSize_AfterInterp)); 
end

ImgSIM_HR_SeparatedBand_Selected=ImgSIM_HR_SeparatedBand(:,:,:,[1:5,7:10,12:15]); % 6, 11 are the zero-order frequency components that are solved repetitively 
ImgSIM_HR_Sum=abs(sum(ImgSIM_HR_SeparatedBand_Selected,4));
FT_ImgSIM_HR_Sum=fftshift(fftn(ifftshift(ImgSIM_HR_Sum)));
ImgSIM_HR_Sum_Truncated=ImgSIM_HR_Sum(:,:,65:192);

figure
subplot(131);imagesc(squeeze(log(abs(1+FT_ImgSIM_HR_Sum(:,:,128)))));axis equal;colormap(hot)
subplot(132);imagesc(squeeze(log(abs(1+FT_ImgSIM_HR_Sum(:,256,:))))');axis equal;colormap(hot)
subplot(133);imagesc(squeeze(log(abs(1+FT_ImgSIM_HR_Sum(256,:,:))))');axis equal;colormap(hot)

figure
subplot(131);imagesc(squeeze(ImgSIM_HR_Sum_Truncated(:,:,64)));axis equal;colormap(hot)
subplot(132);imagesc(squeeze(ImgSIM_HR_Sum_Truncated(:,256,:))');axis equal;colormap(hot)
subplot(133);imagesc(squeeze(ImgSIM_HR_Sum_Truncated(256,:,:))');axis equal;colormap(hot)


% ImgSIM_HR_SeparatedBand_Selected=ImgSIM_HR_SeparatedBand(:,:,:,[1:5]);
% ImgSIM_HR_Sum=sum(ImgSIM_HR_SeparatedBand_Selected,4);
% FT_ImgSIM_HR_Sum=fftshift(fftn(ifftshift(ImgSIM_HR_Sum)));
% figure
% subplot(131);imagesc(squeeze(log(abs(1+FT_ImgSIM_HR_Sum(:,:,128)))));axis equal;colormap(hot)
% subplot(132);imagesc(squeeze(log(abs(1+FT_ImgSIM_HR_Sum(:,256,:))))');axis equal;colormap(hot)
% subplot(133);imagesc(squeeze(log(abs(1+FT_ImgSIM_HR_Sum(256,:,:))))');axis equal;colormap(hot)

