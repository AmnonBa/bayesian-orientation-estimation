function mu = get_sph_density(theta, phi, sph_coef, L)


real_t_complex = get_view_real_t_complex(L);
sph_coef=real_t_complex*sph_coef;


vec_idx = @(p, u) 2*p^2+p+u+1;
mu=0;
for l=0:L
    for m=-2*l:2*l
        Ylm=harmonicY(2*l,m,theta,phi);
        mu=mu+sph_coef(vec_idx(l,m))*Ylm;
    end
end
mu=real(mu);

end