function f = get_vMF_density(mus,w,x,k)
n = size(mus,2);
f = 0;
for i=1:n
    mu = mus(:,i);
    f = f+w(i)*exp(k*sum(mu.*x))*k/2/pi/(exp(k)-exp(-k));
end

if size(f,2)>size(f,1)
f=f';
end
end