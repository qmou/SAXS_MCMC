function [engyf] = mc_move(nx, xsed1, xsed2)

	global mtx_xtal occ_bkb
	global xtal_pos xtal_sidex xtal_sidey xtal_bkbexcl scl
	global n_acpt a_slot kk kr a_max
	%move xtal 1
        engyf = 10000;
        divx = 16;
        divy = 16;
    
	pos1x = xtal_pos(xsed1,1);
	pos1y = xtal_pos(xsed1,2);
	dx = ceil((2*rand(1)-1)/2*nx/divx);	%dx > dy
	dy = ceil((2*rand(1)-1)/2*nx/divy);	%rand between -0.5 +0.5
    
	pos1x_new = modnozero(double(pos1x+dx),nx);
	pos1y_new = modnozero(double(pos1y+dy),nx);

	lmtx = xtal_sidex(xsed1); %sidex = half of actual side length
	lmty = xtal_sidey(xsed1);
    
	pos2x = xtal_pos(xsed2,1);
	pos2y = xtal_pos(xsed2,2);
	dx = ceil((2*rand-1)/2*nx/divx);	%dx > dy
	dy = ceil((2*rand-1)/2*nx/divy);	%rand between -0.5 +0.5

	pos2x_new = modnozero(double(pos2x+dx),nx);
	pos2y_new = modnozero(double(pos2y+dy),nx);

	lmtxp = xtal_sidex(xsed2); %sidex = half of actual side length
	lmtyp = xtal_sidey(xsed2);
    
    r = sqrt(double((pos1x_new-pos2x_new)^2+(pos1y_new-pos2y_new)^2));

	if (occ_bkb(pos1x_new, pos1y_new) == 0) && (occ_bkb(pos2x_new, pos2y_new) == 0) .../
		&& (occ_bkb(modnozero(pos1x_new-lmtx,nx),modnozero(pos1y_new-lmty,nx)) == 0) .../
        && (occ_bkb(modnozero(pos1x_new+lmtx,nx),modnozero(pos1y_new+lmty,nx)) == 0) .../
        && (occ_bkb(modnozero(pos1x_new-lmtx,nx),modnozero(pos1y_new+lmty,nx)) == 0)  .../
        && (occ_bkb(modnozero(pos1x_new+lmtx,nx),modnozero(pos1y_new-lmty,nx)) == 0) .../
        && (occ_bkb(modnozero(pos2x_new+lmtxp,nx),modnozero(pos2y_new-lmtyp,nx)) == 0) .../
        && (occ_bkb(modnozero(pos2x_new+lmtxp,nx),modnozero(pos2y_new+lmtyp,nx)) == 0) .../
        && (occ_bkb(modnozero(pos2x_new-lmtxp,nx),modnozero(pos2y_new-lmtyp,nx)) == 0) .../
        && (occ_bkb(modnozero(pos2x_new-lmtxp,nx),modnozero(pos2y_new+lmtyp,nx)) == 0) .../
        && (r > 2.5*max(lmtx,lmtxp))

        ttxtal_mtx = mtx_xtal;	%temp matrix to hold structure
        
        for k = -lmty:lmty      %paint crystal1
            kyplusk=modnozero(pos1y_new+k,nx);
            for j = -lmtx:lmtx
                    jxplusj = modnozero(pos1x_new+j,nx);
                      	ttxtal_mtx(jxplusj,kyplusk)=0;
            end
    	end

    	for j = -lmtx:lmtx      %unpaint crystal1
            jxplusj = modnozero(pos1x+j,nx);
            for k = -lmty:lmty
                    kyplusk=modnozero(pos1y+k,nx);
                      	ttxtal_mtx(jxplusj,kyplusk)=1;
            end
    	end

    	for j = -lmtxp:lmtxp      %paint crystal2
            jxplusj = modnozero(pos2x_new+j,nx);
            for k = -lmtyp:lmtyp
                    kyplusk=modnozero(pos2y_new+k,nx);
                      	ttxtal_mtx(jxplusj,kyplusk)=0;
            end
    	end

    	for j = -lmtxp:lmtxp      %unpaint crystal2
            jxplusj = modnozero(pos2x+j,nx);
            for k = -lmtyp:lmtyp
                    kyplusk=modnozero(pos2y+k,nx);
                      	ttxtal_mtx(jxplusj,kyplusk)=1;
            end
    	end
	
	% MCMC
    	eng_i = post_p_test(mtx_xtal, scl, 1.85);
    	eng_f = post_p_test(ttxtal_mtx, scl, 1.85);
	
	rto = eng_f/eng_i; %Q(x->x') is symetrical
    
    
	acpt = min(1,1/rto);   %acpt < 1 if the move increase energy
	a_slot(modnozero(kk,200)) = rto;

    	allow = 0;
	
	if rto < 1
		allow = 1;
	end

	if acpt < 1 && rand < 0.2*gaus(acpt, 1, a_max-1 )
		allow = 1;
		kr = kr + 1;
	end

    	if allow > 0
		n_acpt = n_acpt + 1;
    		mtx_xtal = ttxtal_mtx;

    		exlcu1x = xtal_bkbexcl(xsed1,1);
    		exlcu1y = xtal_bkbexcl(xsed1,2);
    		exlcu2x = xtal_bkbexcl(xsed2,1);
    		exlcu2y = xtal_bkbexcl(xsed2,2);

    		for j = -exlcu1x:exlcu1x      %unpaint crystal1 exclu
            	jxplusj = modnozero(pos1x+j,nx);
            	for k = -exlcu1y:exlcu1y
                    kyplusk=modnozero(pos1y+k,nx);
                            occ_bkb(jxplusj,kyplusk)=occ_bkb(jxplusj,kyplusk)-1;
            	end
    		end

    		for j = -exlcu1x:exlcu1x      %paint crystal1 exclu
            	jxplusj = modnozero(pos1x_new+j,nx);
            	for k = -exlcu1y:exlcu1y
                    kyplusk=modnozero(pos1y_new+k,nx);
                      	occ_bkb(jxplusj,kyplusk)=occ_bkb(jxplusj,kyplusk)+1;
            	end
    		end

    		for j = -exlcu2x:exlcu2x      %unpaint crystal1 exclu
            	jxplusj = modnozero(pos2x+j,nx);
            	for k = -exlcu2y:exlcu2y
                    kyplusk=modnozero(pos2y+k,nx);
                            occ_bkb(jxplusj,kyplusk)=occ_bkb(jxplusj,kyplusk)-1;
            	end
    		end

    		for j = -exlcu2x:exlcu2x      %paint crystal1 exclu
            	jxplusj = modnozero(pos2x_new+j,nx);
            	for k = -exlcu2y:exlcu2y
                    kyplusk=modnozero(pos2y_new+k,nx);
                      	occ_bkb(jxplusj,kyplusk)=occ_bkb(jxplusj,kyplusk)+1;
            	end
    		end
    		engyf = eng_f;
            xtal_pos(xsed1,1) = pos1x_new;
            xtal_pos(xsed1,2) = pos1y_new;
            xtal_pos(xsed2,1) = pos2x_new;
            xtal_pos(xsed2,2) = pos2y_new;
        end
    end

end
