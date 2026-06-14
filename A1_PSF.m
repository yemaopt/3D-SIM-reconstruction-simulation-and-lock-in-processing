% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab, 2020.08.20
% Department of Biomedical Engineering, Johns Hopkins University

clear
clc
close all

PixelSize=0.08; %um
z=PixelSize.*(-128:1:127); %um
XYFOV=PixelSize.*128; %um
L=length(z);
PSF=zeros(L+1,L+1,L);

%% Excitation spot with (left) circular polarization
Intensity.Name='Uniform';
PhasePlate.Name='Null';
Polarization='Left circular'; %has tested that left and right circular polarization lead to the same result
for index=1:length(z)
    z(index)
    [ I,X,Y ] = FocusedPatternXY_CZT_RealUnit(Intensity,PhasePlate,Polarization,z(index),XYFOV,PixelSize);
    PSF(:,:,index)=I;
end

PSF=PSF(1:L,1:L,1:L);
PSF=PSF./max(PSF(:)).*65535;
PSF=uint16(PSF);

%%
% Display the results
figure(101);subplot(121);imagesc(PSF(:,:,L/2));axis equal;axis off;subplot(122);imagesc(squeeze(PSF(:,L/2,:))');axis equal;axis off;
% save the results
save('PSF_256Pixels.mat','PSF');
