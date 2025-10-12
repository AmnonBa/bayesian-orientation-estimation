function Psi = precompute_rot_density(SO3_rule, P)

[real_t_complex,~] = get_view_real_t_complex(P);
vec_idx = @(p, u) 2*p^2+p+u+1;

nB = size(real_t_complex,1);
nq = size(SO3_rule,1);
Psi = zeros(nq, nB);

for i=1:nq

    alpha = SO3_rule(i,1);
    beta = SO3_rule(i,2);
    gamma = SO3_rule(i,3);

    for p=0:P
        D_p = wignerD(2*p, alpha, beta, gamma); 


        for u=-(2*p):(2*p)
            Psi(i,vec_idx(p,u)) = D_p(u+2*p+1,0+2*p+1);
        end
    end
end

Psi=real(Psi*real_t_complex);


% [real_t_complex,~] = get_view_real_t_complex(P);
% vec_idx = @(p, u) p^2+p+u+1;
% 
% nB = size(real_t_complex,1);
% nq = size(SO3_rule,1);
% Psi = zeros(nq, nB);
% 
% for i=1:nq
% 
%     alpha = SO3_rule(i,1);
%     beta = SO3_rule(i,2);
%     gamma = SO3_rule(i,3);
% 
%     for p=0:P
%         Wig_p = wignerD(p, alpha, beta, gamma); 
% 
% 
%         for u=-p:p
%             Psi(i,vec_idx(p,u)) = Wig_p(u+p+1,0+p+1);
%         end
%     end
% end
% 
% Psi=real(Psi*real_t_complex);

end 
