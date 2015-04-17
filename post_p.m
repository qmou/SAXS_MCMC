function [resd] =  post_p(asize, cbk, crod, cxtal)

    global nx mtx_xtal mtx_rod
	mtx_rod_d = double(mtx_rod);
	mtx_xtal_d = double(mtx_xtal);
	rho=zeros(nx,nx);
	bk = mtx_rod_d.*mtx_xtal_d;
	rho = cbk*bk + crod*((1-mtx_rod_d).*mtx_xtal_d) + cxtal*(1-mtx_xtal_d);  %the resulting density map in 3D

    [In, Iq] = SAXSrodcylformfac(rho);
    
    if nx == 256
        rdc = 1;
    elseif nx == 1024
        rdc = 3;
    elseif nx == 2048
        rdc = 3;
    elseif nx == 4096
        rdc = 5;
    end

    temp = smoothlog(In,rdc);       
    dq = 2^rdc*2*pi/(nx*asize); %smooth twice, X4
    l = length(temp);
    q = 1:length(temp);
    q = q * dq;

    load 'silk_profile'

	q_spl = 0.16:0.01:1.1;
	aa_spl_y = spline(aa_x,aa_y,q_spl);
	aa_spln_y = aa_spl_y/aa_spl_y(1);

	temp_spl = spline(q,temp,q_spl);
	temp_spln = temp_spl/temp_spl(1);

%     figure
%     loglog(q_spl, aa_spln_y, q_spl, temp_spln);
% 
%     x = nx/2;
%     figure
%     dx = 128;
%     pcolor(rho(x-2*dx:x+2*dx,x-2*dx:x+2*dx));
%     figure
%     pcolor(log(Iq(x-2*dx:x+2*dx,x-dx/2:x+dx/2)));
%     figure
%     loglog(q,temp);

    resd = sum((temp_spln - aa_spln_y).^2);

end