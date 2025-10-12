function Rots = sample_orientations(N,viewing_density,M)


Rots = zeros(3,3,N);

for i=1:N

    while true 
        alpha = rand*2*pi;
        beta = acos(1-2*rand);

        m = viewing_density(beta,alpha);
        if rand<(m/M)
            break
        end
    end

    gamma = rand*2*pi; 
    Rots(:,:,i) = Rz(alpha)*Ry(beta)*Rz(gamma);
end


end