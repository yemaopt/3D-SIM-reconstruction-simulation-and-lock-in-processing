% Written by Ye Ma, PhD student in Prof. Taekjip Ha's lab, 2020.08.20
% Department of Biomedical Engineering, Johns Hopkins University


% clear
% clc
% close all

N=256*3;

num_tubule=1000;
num_point_seed=20;
mol_mu=10;
mol_sigma=3;

Sample_Orig=zeros(N,N,N);

for jj=1:num_tubule
    
    m = log((mol_mu^2)/sqrt(mol_sigma^2+mol_mu^2));
    s = sqrt(log(mol_sigma^2/(mol_mu^2)+1));
    
    pos = zeros(num_point_seed,3);    
    pos(:,1) = rand(num_point_seed,1);
    pos(:,2) = rand(num_point_seed,1);
    pos(:,3) = rand(num_point_seed,1);

    t=linspace(min(min(pos)),max(max(pos)),20);
    x1=spline(t,pos(:,1));
    y1=spline(t,pos(:,2));
    z1=spline(t,pos(:,3));
    
    t1=min(min(pos)):0.0001:max(max(pos));
    x3=ppval(x1,t1);
    y3=ppval(y1,t1);
    z3=ppval(z1,t1);
    
    num_point=length(t1);
    molNum = floor(lognrnd(m,s,num_point,1));

    pos_2 = zeros(num_point,3);
    pos_2(:,1) = min(N,max(1,floor(1+(N-1).*x3)));
    pos_2(:,2) = min(N,max(1,floor(1+(N-1).*y3)));
    pos_2(:,3) = min(N,max(1,floor(1+(N-1).*z3)));

    for ii=1:num_point
        Sample_Orig(pos_2(ii,1),pos_2(ii,2),pos_2(ii,3))=molNum(ii);
    end 
    
%     figure
%     scatter3(pos_2(:,1),pos_2(:,2),pos_2(:,3));
%     grid on
%     xlabel('x/m')
%     ylabel('y/m')
%     zlabel('z/m')
    
end

ll=(-128:127)+N/2;
Sample=Sample_Orig(ll,ll,ll);

figure
imagesc(Sample(:,:,128))

save Sample1000Tubule_3D_256Pixels.mat Sample