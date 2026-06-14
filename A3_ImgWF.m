% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab, 2020.08.20
% Department of Biomedical Engineering, Johns Hopkins University

close all
%clear

load Sample1000Tubule_3D_256Pixels.mat
load PSF_256Pixels.mat
fSample=fftshift(fftn(ifftshift(Sample)));
%%
fPSF=fftshift(fftn(ifftshift(PSF)));
ImgWF=fftshift(ifftn(ifftshift(fSample.*fPSF)));

ImgWF_XY=squeeze(ImgWF(:,:,128));
ImgWF_XZ=squeeze(ImgWF(128,:,:))';
ImgWF_XY=ImgWF_XY./max(max(ImgWF_XY));
ImgWF_XZ=ImgWF_XZ./max(max(ImgWF_XZ));

figure
imagesc(ImgWF_XY);colormap(hot);axis equal;axis off
figure
imagesc(ImgWF_XZ);colormap(hot);axis equal;axis off

fImgWF=fftshift(fftn(ifftshift(ImgWF)));
figure
imagesc(log(abs(squeeze(fImgWF(:,:,128)))+1));colormap(hot);axis equal;axis off
figure
imagesc(log(abs(squeeze(fImgWF(:,128,:)))+1));colormap(hot);axis equal;axis off

