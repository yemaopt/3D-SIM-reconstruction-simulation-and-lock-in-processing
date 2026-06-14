% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab, 2020.08.20
% Department of Biomedical Engineering, Johns Hopkins University

%close all
clc
clear


%refractive index is considered to be 1 here

L=256;
IncidentAngle = pi/3;
ccdPxlSize = 0.08; % in um
emWavelength = 0.488; % in um

[xxx,yyy,zzz] = meshgrid(-128:1:127);
xxx=xxx.*ccdPxlSize;
yyy=yyy.*ccdPxlSize;
zzz=zzz.*ccdPxlSize;

k = 2*pi/emWavelength;
kr = k*sin(IncidentAngle);
kz = k*cos(IncidentAngle);
A1=1;A2=1;A3=1;

Pattern3DSIM=zeros(L,L,L,15);

ii=0;
for theta = [0, pi/3, 2*pi/3]
    for phi1 = 2*pi/5.*(0:1:4)
        ii = ii+1;
        kx = kr*cos(theta);
        ky = kr*sin(theta);
        phi2=phi1/2;
        
        E1=A1.*exp(1j.*(kx.*xxx+ky.*yyy+phi1)+1j.*kz.*zzz);
        E2=A2.*exp(1j.*(k.*zzz+phi2));
        E3=A3.*exp(-1j.*(kx.*xxx+ky.*yyy)+1j.*kz.*zzz);
        Pattern3DSIM(:,:,:,ii)=abs(E1+E2+E3).^2; 
        %         figure
        %         subplot(131);imagesc(squeeze(Pattern3DSIM(:,:,128,ii)));axis equal
        %         subplot(132);imagesc(squeeze(Pattern3DSIM(:,128,:,ii)));axis equal
        %         subplot(133);imagesc(squeeze(Pattern3DSIM(128,:,:,ii)));axis equal
        %         pause()
    end
end

load Sample1000Tubule_3D_256Pixels.mat
load PSF_256Pixels.mat

temp=zeros(L,L,L);
temp(:,:,65:192)=Sample(:,:,65:192);

ImgSIM=zeros(L,L,L,15);
for z0 = -64:1:63
    Sample1=circshift(temp,z0,3);
    z0
    for ii=1:15
        Object = Sample1.*Pattern3DSIM(:,:,:,ii);
        fObject=fftshift(fftn(ifftshift(Object)));
        fPSF=fftshift(fftn(ifftshift(PSF)));
        Img3D=fftshift(ifftn(ifftshift(fObject.*fPSF)));
        ImgSIM(:,:,z0+129,ii)=Img3D(:,:,128);
    end
end

ImgSIM=uint16(ImgSIM./max(ImgSIM(:)).*65535);

save Img3DSIM.mat ImgSIM

% ii=6;
% Object = temp.*Pattern3DSIM(:,:,:,ii);
% fObject=fftshift(fftn(ifftshift(Object)));
% fPSF=fftshift(fftn(ifftshift(PSF)));
% ImgWF=fftshift(ifftn(ifftshift(fObject.*fPSF)));
% fImgWF=fftshift(fftn(ifftshift(ImgWF)));
% 
% figure;
% subplot(221);imagesc(squeeze(ImgWF(:,:,128)));colormap(hot);axis equal;axis off
% subplot(222);imagesc(squeeze(ImgWF(:,128,:))');colormap(hot);axis equal;axis off
% subplot(223);imagesc(log(abs(squeeze(fImgWF(:,:,128)))+1));colormap(hot);axis equal;axis off
% subplot(224);imagesc(log(abs(squeeze(fImgWF(:,128,:))')+1));colormap(hot);axis equal;axis off
% 
figure;ii=1;
subplot(221);imagesc(squeeze(ImgSIM(:,:,128,ii)));colormap(hot);axis equal;axis off
subplot(222);imagesc(squeeze(ImgSIM(:,128,:,ii))');colormap(hot);axis equal;axis off
fImgSIM=fftshift(fftn(ifftshift(ImgSIM(:,:,:,ii))));
subplot(223);imagesc(log(abs(squeeze(fImgSIM(:,:,128)))+1));colormap(hot);axis equal;axis off
subplot(224);imagesc(log(abs(squeeze(fImgSIM(:,128,:))')+1));colormap(hot);axis equal;axis off
