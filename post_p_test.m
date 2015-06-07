function [resd] =  post_p_test(ttmtx, asize, cxtal)

    global nx q_spl aa_spln_y

    [In, ~] = SAXSrodcylformfac(double(cxtal*(1-ttmtx)));
    
    switch nx
%         case 256
%             rdc = 1;
%         case 1024
%             rdc = 3;
        case 2048
            rdc = 3;
        case 4096
            rdc = 4;
    end

    temp = smoothlog(In,rdc);       
    dq = 2^rdc*2*pi/(nx*asize); %smooth twice, X4
    q = 1:length(temp);
    q = q*dq;

	temp_spl = spline(q,temp,q_spl);
	temp_spln = temp_spl/temp_spl(1);

    resd = sum((temp_spln - aa_spln_y).^2);

end