function [resd] =  plot_density(asize, cbk, crod, cxtal)

    global nx mtx_xtal mtx_rod q_spl aa_spln_y occ_bkb xtal_pos
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

	temp_spl = spline(q,temp,q_spl);
	temp_spln = temp_spl/temp_spl(1);

    figure
    loglog(q_spl, aa_spln_y, q_spl, temp_spln);

    figure
    dx = 512;
    pcolor(rho(nx-dx:nx,nx-dx:nx));
    colormap(bone)
    
    figure
    dx = 512;
    occ = double(occ_bkb);
    pcolor(occ(nx-dx:nx,nx-dx:nx));
    
    
    lgh = nx;
    
    for cnt = 1:4
    rd_mtx = zeros(lgh/2,lgh/2);
        for kk = 1:2:lgh
            k = (kk+1)/2;
            for jj = 1:2:lgh
                j = (jj+1)/2;
                temp = rho(kk,jj) + rho(kk,jj+1) + rho(kk+1,jj) + rho(kk+1,jj+1);
                temp = temp / 4;
                rd_mtx(k,j) = temp;
            end
        end
    clear rho
    rho = rd_mtx;
    lgh = lgh / 2;
    end
    
    figure
    pcolor(rd_mtx);
    lmt = 170;
    rho_0 = sum(sum(1-mtx_xtal))/(nx^2);
    pdf_aa = rdf(xtal_pos);
    l_aa = 1:length(pdf_aa);
    l_aa = l_aa*asize;
    for i = 1:length(pdf_aa)
        pdf_aa(i) = pdf_aa(i)/(4*pi*(l_aa(i)^2)*rho_0);
    end
    figure
    plot(l_aa(1:lmt),pdf_aa(1:lmt),'LineWidth',1.2);
    resd = sum((temp_spln - aa_spln_y).^2);

end