clear all
global nx mtx_xtal mtx_rod occ_xtal occ_bkb egyf
global xtal_pos xtal_sidex xtal_sidey xtal_bkbexcl scl
global n_acpt
global a_slot kk kr a_max
nx=2048;
rg_size = 128;

scl = 0.305;

xtal_pos = int16([]);
xtal_sidex = int16([]);
xtal_sidey = int16([]);
xtal_bkbexcl = int16([]);
%% crystal parameter
disp('---> define crystal phase');
xtal_life = 1000000;   %squre crystal count
if xtal_life > 0
    xtal_minlength = 6;
    xtal_maxlength = 8;
    xtal_minlength = xtal_minlength/2;
    xtal_maxlength = xtal_maxlength/2;
    xtal_range = xtal_maxlength - xtal_minlength;
    xtal_lmt = 0.25;    %crystal fraction
    xtal_asqmin = 1;
    xtal_asqmax = 1.1;
    axtal_range = xtal_asqmax-xtal_asqmin;
%     xtal_exclu_inter = 3.2;
    xtal_exclu_inter = 4.0;
    xtal_eocc_bkb = 1.6; %was 2.2
end

occ_xtal = int8(zeros(nx,nx));   %book keeping of excluded region, hard shell replusive model
occ_bkb = int16(zeros(nx,nx));      %exclude region for backbone, more relax than occ_xtal
%%
disp('---> define amorphous phase, long rod configuration');
rod_life = 0;   %long rod number
if rod_life > 0
    rod_radmin = 1;
    rod_radmin = rod_radmin/2;
    rod_radmax = 3;
    rod_radmax = rod_radmax/2;
    rod_range = rod_radmax - rod_radmin;
    rod_aspmin = 8;
    rod_aspmax = 15;
    rod_asprange = rod_aspmax-rod_aspmin;
    rod_lmt = 0.15;
    rod_exclu = 3;
    rod_exclu_y = 0.6;
end

%%
rg_cnt = nx/rg_size;
rg_cnt_half = rg_cnt/2;
%%
disp('---> generating crystal phase');
mtx_rod=int8(ones(nx,nx));  
mtx_xtal = int8(ones(nx,nx));

n_xtal = 0;
xtal_frac = 0;

for cnt = 1:rg_cnt
    if mod(cnt,2) == 1     %skip some number region
        for nc = 1:xtal_life
                if xtal_frac < xtal_lmt
                    rg_shift = rg_size*(cnt-1);    %cal real shift of block from origin
                    jx = ceil(rand*rg_size);      %confine to the near center region
                    ky = ceil(rand*nx);      
                    allowed = 0;
                    if occ_xtal(jx+rg_shift,ky)<1
                        allowed = 1;
                        xtal_pos = cat(1,xtal_pos,int16([jx+rg_shift,ky]));
                    end
                    if allowed == 1
                        side_x = xtal_minlength + rand*xtal_range;
                        n_xtal = n_xtal + 1;
                        side_x = round(side_x);
                        side_y = round(side_x*(xtal_asqmin + rand*axtal_range));
                        xtal_sidex = cat(1,xtal_sidex,side_x);
                        xtal_sidey = cat(1,xtal_sidey,side_y);

                        paint_xtal(jx, ky, side_x, side_y, rg_shift, 'xtal');
                        
                        sqxlulength = round(side_x*xtal_exclu_inter);
                        sqxlulengthy = round(side_y*xtal_exclu_inter);
                        
                        paint_xtal(jx, ky, sqxlulength, sqxlulengthy, rg_shift, 'occxtal');

                        %exclude matrix for backbone structure, smaller exclusion
                        %factor
                        sqxlulength = round(side_x*xtal_eocc_bkb);
                        sqxlulengthy = round(side_y*xtal_eocc_bkb);

                        xtal_bkbexcl = cat(1, xtal_bkbexcl, [sqxlulength,sqxlulengthy]);

                        paint_xtal(jx, ky, sqxlulength, sqxlulengthy, rg_shift, 'occbkb');

                        xtal_frac=1-double(sum(sum(mtx_xtal)))/double(nx^2);
                    end     %if allowed = 1
                end     %end xtal_frac< fsq
        end
    end
end

rg_shift = 1;
%% generate amophous phase
disp('---> generating amorphous phase');
clustmult = 10;
n_rod = 0;
rod_frac = 0;
lmt = 3;
rod_life_rdc = ceil(rod_life/clustmult);
for nc = 1:clustmult
    for n = 1:rod_life_rdc  %otherwise 1:rod_life takes too much memory
        if rod_frac < rod_lmt  %stop if crystallinity is reached
             jx = ceil(rand*nx);  
             ky = ceil(rand*nx);
             allowed = 0;       
             if occ_bkb(jx,ky) < 1 .../
                     && occ_bkb(modnozero(jx-2*lmt,nx),modnozero(ky-lmt,nx)) == 0 .../
                     && occ_bkb(modnozero(jx+2*lmt,nx),modnozero(ky+lmt,nx)) == 0 .../
                     && occ_bkb(modnozero(jx-2*lmt,nx),modnozero(ky+lmt,nx)) == 0 .../
                     && occ_bkb(modnozero(jx+2*lmt,nx),modnozero(ky-lmt,nx)) == 0
                 allowed = 1;
             end %if occ_bkb
             if allowed == 1
                  n_rod = n_rod + 1;
                  side_x = rod_radmin+rand*rod_range; %normal distribution with cut-off radmin
                  nside_x = round(side_x);
                  aspratio = rod_aspmin+rand*rod_asprange;
                  nside_y = round(side_x*aspratio);     % x=z dimension, modify v1, change ratio of Z
                  paint_xtal(jx, ky, nside_x, 2*nside_y, rg_shift, 'rod');
                  
                  nradexcl = round(nside_x*rod_exclu);
                  nradexcly = round(nside_x*(rod_exclu - rod_exclu_y)); %modifier, reduce exclusion
                  paint_xtal(jx, ky, nradexcl, nradexcly, rg_shift, 'rodocc')

                  rod_frac=1-double(sum(sum(mtx_rod)))/double(nx^2);  %occ_bkbupied fraction
              end  %if allowed
        end %if <fcryst
    end  %n clusters
end  %rod_lifemult
%%
% assign contribution of diffraction for different structures

load 'silk_pfclean'
global q_spl aa_spln_y

q_spl = 0.155:0.005:1.2;
aa_spl_y = spline(nc_new_x, nc_new_y, q_spl);
aa_spln_y = aa_spl_y/aa_spl_y(1);

xtal_pop = length(xtal_pos);
egyf = post_p_test(mtx_xtal, scl, 1.85);

kk = 0;
n_acpt = 0;
a_slot = zeros(200,1);
%a_min = 0.96;
a_max = 1.015;
kr = 0;

while egyf > 1e-3
	xsed1 = floor(rand*xtal_pop);
	xsed2 = floor(rand*xtal_pop);
    
    if modnozero(kk,2000) == 2000
        a_max = max(a_slot);
        if a_max < 1.01
            a_max = 1.01;
        end
        if a_max <= 1.0
            a_max = 1.005;
        end
    end
    
	if xsed1 == 0
		xsed1 = floor(rand*xtal_pop/2)+1;
	end
	if xsed2 == 0
		xsed2 = floor(rand*xtal_pop/2)+1;
	end

	tmp = mc_move(nx,xsed1,xsed2);
    
    if tmp < egyf
        egyf = tmp;
    end
    
    if mod(kk,2000) == 0
        disp(sprintf('[Engy]%7.5f, [Int]%d, [a_max]%4.3f, [a_min]%4.3f, [n_acpt]%d, [mc]%d', egyf, kk, a_max, min(a_slot), n_acpt, kr));
    end
    kk = kk + 1;
    if kk > 10e5
      break
    end
end
sprintf('Accept No. %d, MC No. %d', n_acpt, kr)
