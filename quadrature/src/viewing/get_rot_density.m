function f = get_rot_density(Rots, rot_coef, P, is_real_coef, first_entry_remove)

vec_idx = @(p, u) p^2+u+p+1; 

if first_entry_remove==true 
    rot_coef=[1;rot_coef];
end

if is_real_coef==true 
    [real_t_complex,~] = get_view_real_t_complex(P);
    rot_coef=real_t_complex*rot_coef;
end

if numel(size(Rots))==3
    N = size(Rots,3);
    f = zeros(N,1);
    for i=1:N
        Rot = Rots(:,:,i);
        [alpha, beta, gamma] = rot2elr(Rot);
    
        Psi = zeros(numel(rot_coef),1);
        for p=0:P
            D_p = wignerD(p, alpha, beta, gamma); 
        
        
            for u=-p:p
                Psi(vec_idx(p,u)) = D_p(u+p+1,0+p+1);
            end
        end
        
        f(i) = real(sum(Psi.*rot_coef));
    end

else

        Rot = Rots;
        [alpha, beta, gamma] = rot2elr(Rot);
    
        Psi = zeros(numel(rot_coef),1);
        for p=0:P
            D_p = wignerD(p, alpha, beta, gamma); 
        
        
            for u=-p:p
                Psi(vec_idx(p,u)) = D_p(u+p+1,0+p+1);
            end
        end
        
        f = real(sum(Psi.*rot_coef));
end

end