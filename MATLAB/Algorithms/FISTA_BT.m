function [res,stepvals] = FISTA_BT(proj,geo,angles,niter,hyper,eta)
% FISTA is a quadratically converging algorithm that relies on the hyper
% parameter named 'hyper'. This parameter should approximate the largest 
% eigenvalue in the A matrix in the equation Ax-b and Atb. Empirical tests
% show that for, the headphantom object:
%           geo.nVoxel = [64,64,64]'    ,      hyper (approx=) 2.e8
%           geo.nVoxel = [512,512,512]' ,      hyper (approx=) 2.e4
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD.
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri, Reuben Lindroos
%--------------------------------------------------------------------------

lambda = 0.01;
res = zeros(geo.nVoxel','single');
x_rec = res;
L = hyper;
stepvals = []

oldDiff = 0;
for ii = 1:niter
    Lbar = L; 
    if (ii == 1); 
        disp('Estimating maximum value of L.');
        jj = 1;
        while true
            tic;
            lambda0 = lambda/Lbar; 
            zk = res - 1/Lbar*Atb(proj - Ax(res,geo,angles),geo,angles);
            calc_f = proj-Ax(zk,geo,angles);
            F = 0.5*norm(calc_f(:)) + norm(lambda0*res(:),1);
            Q = calc_Q(proj, geo, angles, zk, res, Lbar,lambda0);
            if mod(jj,2) == 0;
                fprintf('L: %e | iteration: %d | Time(s)/iter: %.2f\n',L, jj,toc)
            end
            if F <= Q
                break;
            end
            Lbar = Lbar*eta; 
            L = Lbar;
            jj = jj +1;
        end
    end
    
    bm = 1/L;
    t = 1;
    lambdaforTV = 2* bm* lambda;
    
    fwd_proj = proj - Ax(res,geo,angles, 'ray-voxel');
    newStepVal = lambdaforTV * im3Dnorm(res,'TV','forward') +norm(fwd_proj(:),2);
    
    if ii == 1; oldStepVal = newStepVal;end
    
    newDiff = abs(oldStepVal - newStepVal)/newStepVal;
    if newDiff >= oldDiff
        
        L = L*0.9;
        fprintf('L: %e | ii: %d | old: %f | new: %f \n',L, ii,oldDiff,newDiff);
        bm = 1/L;
        lambdaforTV = 2* bm* lambda;
    elseif newDiff <=0.95*oldDiff
        L = L*1.1
        bm = 1/L;
        lambdaforTV = 2* bm* lambda;
    end
    stepvals = [stepvals;newStepVal];
    oldDiff = newDiff;
    oldStepVal = newStepVal;
    % gradient descent step
    res = res + bm * 2 * Atb(fwd_proj, geo, angles, 'matched');
    x_recold = x_rec;
    x_rec = im3DDenoise(res,'TV',20,1/lambdaforTV);  
    told = t;
    t = ( 1+sqrt(1+4*t*t) ) / 2;
    res= x_rec + (told-1)/t * (x_rec - x_recold);
    if (ii==1);
        expected_time=toc*niter;
        disp('FISTA');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
    
end


%% compute g
function res_g = g(x,lambda) 
        res_g = sum(abs(lambda*x(:)));
end
%% compute Q
function res_Q = calc_Q(proj,geo,angles,x, y, L,lambda) 
    % based on equation 2.5, page 189
    calc_fy = proj - Ax(y,geo,angles);
    diff = x - y;
    back_projection =  Atb(calc_fy,geo,angles);
    res_Q = 0.5*norm(calc_fy(:)) + dot(diff(:), back_projection(:)) ...
    + L/2*norm(diff(:)) + g(x,lambda);
end

% https://github.com/tiepvupsu/FISTA/blob/master/fista_backtracking.m
end