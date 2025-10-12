function out = form_estimated_moments(V, d, Rots, params)

rng(1)

Nrot = size(Rots,3);
batch_size = 1e3;
D = size(V,1);

n_streams = ceil(Nrot/batch_size);
s2 = params.s2;
s3 = params.s3;
r2_max = params.r2_max ;
r3_max = params.r3_max ;
tol2 = params.tol2;
tol3 = params.tol3;
sigma = params.sigma;
nstd = params.nstd;
augsize = params.augsize;

if augsize>0
    N = 2*augsize*Nrot;
else
    N = Nrot;
end



% random sketch matrix 
G2 = randn(d^2, s2);
G3_1 = randn(d^2, s3);
G3_2 = randn(d^2, s3);

% skecth the moments 
M1 = 0;
Y2 = 0;
Y3 = 0;
time_sketch = 0;
for i=1:n_streams
    disp(['sketch the images over the stream No.',num2str(i), '/', num2str(n_streams)])
    rots = Rots(:,:,((i-1)*batch_size+1):min(i*batch_size,Nrot));
    % tic 
    projections = vol_project(V,rots,nstd,i,augsize);
    % toc

   
    % tic 
    projections = projections_preprocessing(projections,d,sigma);
    % toc


    tic 
    [M1, Y2, Y3] = projections_t_sketched_MoMs(projections, M1, Y2, Y3, G2, G3_1, G3_2, N);
    t=toc;
    time_sketch = time_sketch + t; 
end


H = get_preprocessing_map(D,d,sigma);


% debias
disp(['debiasing after sketching'])
tic 


if nstd>0
    B2 = nstd^2*H*(H'*G2);
    B3 = 0;
    for i=1:D^2
        Hi = H(:,i);
        B3 = B3 + M1*((G3_1'*Hi).*G3_2'*Hi).';
        B3 = B3 + Hi*((G3_1'*M1).*G3_2'*Hi).';
        B3 = B3 + Hi*((G3_1'*Hi).*G3_2'*M1).';
    end
    B3 = nstd^2*B3; 
    
    Y2 = Y2-B2;
    Y3 = Y3-B3;
end 

time_debias1=toc;


tic 
% obtain subspaces by SVD
[U2,S2,~] = svd(Y2,'econ'); 
sig2 = diag(S2);
if sum(sig2(1:r2_max).^2)/sum(sig2.^2)<(1-tol2)
    r2=r2_max; 
else
    r2 = find(cumsum(sig2.^2)/sum(sig2.^2)>(1-tol2),1,'first');
end
U2 = U2(:,1:r2);


[U3,S3,~] = svd(Y3,'econ');
sig3 = diag(S3);
if sum(sig3(1:r3_max).^2)/sum(sig3.^2)<=(1-tol3)
    r3=r3_max; 
else
    r3 = find(cumsum(sig3.^2)/sum(sig3.^2)>(1-tol3),1,'first');
end
U3 = U3(:,1:r3);

time_svd = toc; 

U1= U2;
r1= size(U1,2);
% form the first three subspace moments
m1 = U1'*M1; 
m2 = 0;
m3 = 0;
time_form = 0;
for i=1:n_streams
    disp(['forming the moments over the stream No.',num2str(i), '/', num2str(n_streams)])
    rots = Rots(:,:,((i-1)*batch_size+1):min(i*batch_size,Nrot));
    % tic 
    projections = vol_project(V,rots,nstd,i,augsize);
    % toc 

    % tic 
    projections = projections_preprocessing(projections,d,sigma);
    % toc 

    tic
    [m2, m3] = projections_t_subspace_MoMs(projections, m2, m3, U2, U3, N);
    t=toc;
    time_form = time_form+t;
end


% debias
disp(['debiasing the subspace moments'])
tic 

if nstd>0
    b2 = U2'*H;
    b2 = nstd^2*(b2*b2');

    b3 = 0;
    U3_M1 = U3'*M1;
    for i=1:D^2
        Hi = H(:,i);
        U3_Hi = U3'*Hi;
        b3 = b3 + tns_kron3(U3_M1,U3_Hi,U3_Hi);
        b3 = b3 + tns_kron3(U3_Hi,U3_M1,U3_Hi);
        b3 = b3 + tns_kron3(U3_Hi,U3_Hi,U3_M1);
    end
    b3 = nstd^2*b3; 

    m2 = m2-b2;
    m3 = m3-b3;

end 

time_debias2=toc; 

% subspace moments  
out.m1=d*m1;
out.m2=d^2*m2;
out.m3=d^3*m3;

% bases of subspaces
for i=1:r1
    I = reshape(U1(:,i),[d,d]);
    I = centered_fft2(I)/d;
    U1(:,i) = I(:);
end

for i=1:r2
    I = reshape(U2(:,i),[d,d]);
    I = centered_fft2(I)/d;
    U2(:,i) = I(:);
end

for i=1:r3
    I = reshape(U3(:,i),[d,d]);
    I = centered_fft2(I)/d;
    U3(:,i) = I(:);
end

out.U1=U1;
out.U2=U2;
out.U3=U3;

% singular values 
out.sig2=sig2;
out.sig3=sig3;


% skecthed moments 
out.M1=M1;
out.Y2=Y2;
out.Y3=Y3;

% random sketch operators 
out.G2=G2;
out.G3_1=G3_1;
out.G3_2=G3_2;

% timings
out.time_sketch=time_sketch;
out.time_svd=time_svd;
out.time_form=time_form;
out.time_debias1=time_debias1;
out.time_debias2=time_debias2;


end
%
%
%
%
%
function [M1, Y2, Y3] = projections_t_sketched_MoMs(projections, M1, Y2, Y3, G, G1, G2, N)


for k=1:numel(projections)

    IR = projections{k};
    IR = IR(:);

    M1 = M1 + IR/N;
    Y2 = Y2 + IR*(IR'*G)/N;
    Y3 = Y3 + IR*((G1'*IR).*(G2'*IR)).'/N;

end
    
    
end
%
%
%
%
%
function [m2, m3] = projections_t_subspace_MoMs(projections, m2, m3, U2, U3, N)


for k=1:numel(projections)

    IR = projections{k};
    IR = IR(:);

    I2 = U2'*IR;
    I3 = U3'*IR;

    m2 = m2 + (I2*I2')/N;
    m3 = m3 + tns_kron(tns_kron(I3,I3),I3)/N;

end

end 
