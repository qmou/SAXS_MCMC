function paint_xtal(jx, ky, side_x, side_y, rg_shift, tpy)

	global nx mtx_xtal mtx_rod occ_xtal occ_bkb

	if strcmp(tpy, 'xtal')
		for j = -side_x:side_x      %paint crystal
            jxplusj = jx+j;
            jxplusj = modnozero(jxplusj+rg_shift,nx);
            for k = -side_y:side_y
                    kyplusk=modnozero(ky+k,nx);
                      	mtx_xtal(jxplusj,kyplusk)=0;
            end
    	end

 	elseif strcmp(tpy, 'rod')
 		for j = 0:side_y      %paint rod, inversed order to maintain rod along vertical
            jxplusj=modnozero(jx+j,nx);
            for k = 0:side_x
                    kyplusk=modnozero(ky+k,nx);
                      	mtx_rod(jxplusj,kyplusk)=0;
            end
    	end

    elseif strcmp(tpy, 'occxtal')
    	for j = -side_x:side_x      %paint crystal
            jxplusj = jx+j;
            jxplusj = modnozero(jxplusj+rg_shift,nx);
            for k = -side_y:side_y
                    kyplusk=modnozero(ky+k,nx);
                      	occ_xtal(jxplusj,kyplusk)=1;
            end
    	end

    elseif strcmp(tpy, 'occbkb')
    	for j = -side_x:side_x      %paint crystal
            jxplusj = jx+j;
            jxplusj = modnozero(jxplusj+rg_shift,nx);
            for k = -side_y:side_y
                    kyplusk=modnozero(ky+k,nx);
                      	occ_bkb(jxplusj,kyplusk)=occ_bkb(jxplusj,kyplusk)+1;
            end
    	end

    elseif strcmp(tpy, 'rodocc')
    	for j = -side_x:round(1.5*side_x)      %paint crystal
            jxplusj=modnozero(jx+j,nx);
            for k = -side_y:round(1.5*side_y)
                    kyplusk=modnozero(ky+k,nx);
                      	occ_bkb(jxplusj,kyplusk)=1;
            end
    	end
    end

end