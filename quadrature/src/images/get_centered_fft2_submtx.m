function F2 = get_centered_fft2_submtx(n,row_id,col_id)

if mod(n,2)==0
    k = ((-n/2):(n/2-1))/n;  
else
    k = ((-(n-1)/2):((n-1)/2))/n; 
end
[k1,k2] = meshgrid(k);
k1 = k1(:);
k2 = k2(:);

if numel(row_id)>0
    k1 = k1(row_id);
    k2 = k2(row_id);
end

if mod(n,2)==0
    x = 2*pi*((-n/2):(n/2-1)); 
else
    x = 2*pi*((-(n-1)/2):((n-1)/2)); 
end
[x1,x2] = meshgrid(x);
x1 = x1(:);
x2 = x2(:);

if numel(col_id)>0
    x1 = x1(col_id);
    x2 = x2(col_id);
end

F2 = exp(-1j*(k1*x1'+k2*x2'));

end