function H = get_preprocessing_map(D,d,sigma)

 
filt_mtx = get_filt_mtx(d,sigma);
filt_mtx = sparse(diag(filt_mtx(:)));



ds_start=floor(D/2)-floor(d/2);
ds_idx=(ds_start+1):(ds_start+d);
[ind1,ind2] = ndgrid(ds_idx,ds_idx);
ind1 = ind1(:);
ind2 = ind2(:);
ind = sub2ind([D,D],ind1,ind2);



F = get_centered_fft2_submtx(D,ind,[]);



circ_mask = get_circ_mask(D);
circ_mask = circ_mask(:);
for i=1:d^2
    F(i,:) = F(i,:).'.*circ_mask;
end



H = (d^2/D^2)*filt_mtx*F;


f =  get_centered_fft2_submtx(d,[],[]);
H = f'*H/d^2;


end