function circ_mask = get_circ_mask(n)
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
end