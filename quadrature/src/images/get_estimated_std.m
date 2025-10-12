function std = get_estimated_std(V, SNR)
N=1000;
Rots=zeros(3,3,N);
for i=1:N
    Rots(:,:,i)=get_rand_Rot; 
end
D=size(V,1);
projections=vol_project(V, Rots,0,1,0);
std=0;
for i=1:N
    E=randn(D,D);
    I=projections{i};
    std=std+(sqrt(norm(I,'fro')^2/(SNR*norm(E,'fro')^2)))/N;
end
end