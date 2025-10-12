function b_pu = vMF_t_rot_coef(mus,w,k,P)

if k==0
    func = @(alpha, beta, gamma) 1;
else
    func = @(alpha, beta, gamma) 4*pi*get_vMF_density(mus,w,[sin(beta)*cos(alpha);sin(beta)*sin(alpha);cos(beta)],k);
end

B_puv = WignerD_transform(func, 2*P);
B_puv = B_puv/B_puv(1);

b_pu = B_puv(get_sphere_part(2*P));

idx = 0;
idx_set = [];
for p=0:2*P
    for u=-p:p
        idx = idx+1;
        if mod(p,2)==0
            idx_set = [idx_set, idx];
        end
    end
end
b_pu = b_pu(idx_set);

[~, complex_t_real] = get_view_real_t_complex(P);

b_pu = complex_t_real*b_pu;
b_pu = b_pu(2:end);



% if k==0
%     func = @(alpha, beta, gamma) 1;
% else
%     func = @(alpha, beta, gamma) 4*pi*get_vMF_density(mus,w,[sin(beta)*cos(alpha);sin(beta)*sin(alpha);cos(beta)],k);
% end
% 
% B_puv = WignerD_transform(func, P);
% B_puv = B_puv/B_puv(1);
% 
% b_pu = B_puv(get_sphere_part(P));
% 
% 
% 
% [~, complex_t_real] = get_view_real_t_complex(P);
% if is_real_coef==true
%     b_pu = complex_t_real*b_pu;
% end
% 
% if remove_first_entry==true
%     b_pu = b_pu(2:end);
% end

end




function B_puv = WignerD_transform(func, P)

alpha = pi/(P+1)*(0:(2*P+1));
gamma = pi/(P+1)*(0:(2*P+1));

h_alpha = alpha(2)-alpha(1);
w1 = ones([1, 2*(P+1)])*h_alpha;
h_gamma = gamma(2)-gamma(1);
w3 = ones([1, 2*(P+1)])*h_gamma;

[x,w2]=lgwt(2*(P+1),-1,1);

x = flip(x);
w2 = flip(w2);

B_puv = [];


for l=0:P

    temp = zeros([2*l+1, 2*l+1]);
    
    for j1=1:2*(P+1)
        for k=1:2*(P+1)
            for j2=1:2*(P+1)

                [D_l, ~, ~] = wignerD(l,alpha(j1), acos(x(k)), gamma(j2));
                temp = temp + (2*l+1)*(w1(j1)*w2(k)*w3(j2))*func(alpha(j1), acos(x(k)), gamma(j2))*conj(D_l)/(8*pi^2);

            end
        end
    end

    temp = temp.';
    B_puv = [B_puv; temp(:)];

end

end


function idx = get_sphere_part(P)
idx = [];
count = 1;
for p=0:P
    for u=-p:p
        for v=-p:p
            if v==0
                idx = [idx; count];
            end
            count = count+1;
        end
    end
end
end


