function projections = vol_project(V,Rots,std,seed,augsize)

warning('off')
rng(seed)
n=size(V,1);
accuracy=1e-8;
N = size(Rots,3);


% frequency grids 
if mod(n,2)==0
    k = ((-n/2):(n/2-1))/n;  
else
    k = ((-(n-1)/2):((n-1)/2))/n; 
end
[kx,ky] = meshgrid(k);
kx = kx(:); ky = ky(:); 
Rad = (k(end)-k(1))/2;
circ_mask = sqrt(kx.^2+ky.^2)<=Rad;
circ_mask = reshape(circ_mask,[n,n]);


% rotate grids 
rotated_grids = zeros(3,n^2,N);
for i=1:N
    R = Rots(:,:,i);
    rotated_grids(:,:,i) = R(:,1)*kx'+R(:,2)*ky';
end
s = rotated_grids(1,:,:); 
t = rotated_grids(2,:,:); 
u = rotated_grids(3,:,:); 

S=2*pi*s(:);
T=2*pi*t(:);
U=2*pi*u(:);


IF_rot = finufft3d2(S,T,U,-1,accuracy,pagetranspose(V));
IF_rot = reshape(IF_rot, [n^2, N]);


projections = cell(1,N);
for j=1:N
    I= fftshift(ifft2(ifftshift(reshape(IF_rot(:,j),[n,n]))));
    projections{j} = circ_mask.*(I+std*randn(n,n));
end 


if augsize>0

    angles = 2*pi*(1:(augsize-1))/augsize;
    kx_rotated = zeros(n^2,augsize-1);
    ky_rotated = zeros(n^2,augsize-1);
    for i=1:(augsize-1)
        angle = angles(i);
        kx_rotated(:,i) = 2*pi*(kx*cos(angle)-ky*sin(angle));
        ky_rotated(:,i) = 2*pi*(kx*sin(angle)+ky*cos(angle));
    end
    kx_rotated = kx_rotated(:);
    ky_rotated = ky_rotated(:);

    projections_augmented = cell(1,2*N*augsize);
    for j=1:N
        I = projections{j};
        IF_rot = finufft2d2(kx_rotated,ky_rotated,-1,1e-12,I');
        IF_rot = reshape(IF_rot,[n^2,augsize-1]);
        for i=1:augsize
            if i<augsize
                Ia=fftshift(...
                    ifft2(ifftshift(reshape(IF_rot(:,i),[n,n]))));
                projections_augmented{(j-1)*augsize+i}=circ_mask.*Ia;
            else
                projections_augmented{(j-1)*augsize+i}=I;
            end
        end
    end

    kx_ref = -2*pi*kx;
    ky_ref = 2*pi*ky;
    for j=(N*augsize+1):(2*N*augsize)
        I = projections_augmented{j-N*augsize};
        IF_rot = finufft2d2(kx_ref,ky_ref,-1,1e-12,I');
        I = fftshift(ifft2(ifftshift(reshape(IF_rot,[n,n]))));
        projections_augmented{j}=I;
    end

    projections = projections_augmented;
    % N = numel(projections);
    % 
    % 
    % projections_ref = cell(1,N);
    % kx_ref = -2*pi*kx;
    % ky_ref = 2*pi*ky;
    % 
    % for j=1:N
    %     I = projections{j};
    %     IF_rot = finufft2d2(kx_ref,ky_ref,-1,1e-12,I');
    %     I = fftshift(ifft2(ifftshift(reshape(IF_rot,[n,n]))));
    %     projections_ref{j}=I;
    % end
    % 
    % projections_tol = cell(1,2*N);
    % for j=1:2*N
    %     if j<=N
    %         projections_tol{j} = projections{j};
    %     else
    %         projections_tol{j} =  projections_ref{j-N};
    %     end
    % end
    % 
    % projections = projections_tol;
    % size(projections_tol)
    
end


end



    
    