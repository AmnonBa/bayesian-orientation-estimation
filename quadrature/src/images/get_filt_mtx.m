function filt_mtx = get_filt_mtx(n,sigma)

if sigma==1

    
    filt_mtx = 1;


else

    if mod(n,2)==0
        k = 2*pi*((-n/2):(n/2-1))/n;  
    else
        k = 2*pi*((-(n-1)/2):((n-1)/2))/n; 
    end
    [kx,ky] = meshgrid(k);
    
    filt_mtx = exp(-sigma^2*(kx.^2+ky.^2)/2);

end

end