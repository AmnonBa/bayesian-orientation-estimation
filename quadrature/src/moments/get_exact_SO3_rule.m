function SO3_rule = get_exact_SO3_rule(L)


nx = L+1;
x = linspace(0,2*pi,nx+1);
x = x(1:nx)';
wx = ones(1,nx)'/nx;


ny = ceil((L+1)/2);
[y, wy] = lgwt(ny,-1,1);
wy = wy/2;

z = x; 
wz = wx; 

[alpha,beta,gamma] = ndgrid(x,acos(y),z);
alpha = alpha(:);
beta = beta(:);
gamma = gamma(:);



[wx,wy,wz] = ndgrid(wx,wy,wz);
wx = wx(:);
wy = wy(:);
wz = wz(:);
w = wx.*wy.*wz;


SO3_rule = zeros(numel(w),4);


SO3_rule(:,1) = alpha;
SO3_rule(:,2) = beta;
SO3_rule(:,3) = gamma;
SO3_rule(:,4) = w;


end 