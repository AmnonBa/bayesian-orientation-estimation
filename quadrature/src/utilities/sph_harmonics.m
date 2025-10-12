function Yl = sph_harmonics(l, th, phi)
n = numel(th);
Yl = zeros(n,2*l+1);
for m=-l:l
    Yl(:,m+l+1)=harmonicY(l,m,th,phi);
end

end