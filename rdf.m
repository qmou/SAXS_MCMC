function [pdf] = rdf(xtal_pos)

n = length(xtal_pos);

nmax = int16(ceil(sqrt(2*(2048-1)^2)));

pdf = zeros(nmax,1);

for i = 1:(n-1)
    for j = i+1 : n
        
        r_ix = double(xtal_pos(i,1));
        r_iy = double(xtal_pos(i,2));
        r_jx = double(xtal_pos(j,1));
        r_jy = double(xtal_pos(j,2));
        
        r_ij = sqrt((r_ix-r_jx)^2+(r_iy-r_jy)^2);
        
        if r_ij < 2
            disp(sprintf('r=(%d, %d), id=(%d, %d)', r_ix, r_iy, i, j));
        end
        irabs = round(r_ij);
        ishar = r_ij - irabs;
        shar = abs(ishar);
        isignshar = sign(ishar);
        
        pdf(irabs+1) = pdf(irabs+1) + (1 - shar);
        pdf(irabs+isignshar+1) = pdf(irabs+isignshar+1) + shar;
        
    end
end

end