function vol_coef = reflect_vol_coef(vol_coef, L, S, is_real_coef)
if is_real_coef==true
    [real_t_complex,~] = get_vol_real_t_complex(L,S);
    vol_coef = real_t_complex*vol_coef;
end 

idx=1;
for l=0:L
    for s=1:S(l+1)
        for m=-l:l
            vol_coef(idx)=vol_coef(idx)*((-1)^(l-m));
            idx=idx+1;
        end
    end
end

if is_real_coef==true
    vol_coef=real(real_t_complex\vol_coef);
end 

end 