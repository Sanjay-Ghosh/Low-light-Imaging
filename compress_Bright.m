function g = compress_Bright(src, coeff, sigmaS, K, omega)
% create spatial filter
w  = round(6*sigmaS); if (mod(w,2) == 0); w  = w+1; end
filt     = fspecial('gaussian', [w w], sigmaS);
numer=coeff(1).*imfilter(src, filt, 'symmetric');
denom=coeff(1).*ones(size(src));
% size(coeff)

for k = 2:K
	Iphase = omega*(k-1)*src;
	Ic = cos(Iphase);
	Is = sin(Iphase);
	CIc  = imfilter(Ic, filt, 'symmetric');
	CIs  = imfilter(Is, filt, 'symmetric');
	CIcp = imfilter(Ic.*src, filt, 'symmetric');
	CIsp = imfilter(Is.*src, filt, 'symmetric'); 
  %  (2*k-2)
    numer = numer + coeff(2*k-2)*(Ic.*CIcp+Is.*CIsp);
 	denom = denom + coeff(2*k-2)*(Ic.*CIc +Is.*CIs );

    
% 	Ic = cos(Iphase);
% 	Is = sin(Iphase);
% 	CIc  = imfilter(Ic, filt, 'symmetric');
% 	CIs  = imfilter(Is, filt, 'symmetric');
% 	CIcp = imfilter(Ic.*src, filt, 'symmetric');
% 	CIsp = imfilter(Is.*src, filt, 'symmetric'); 


%     numer = numer + coeff(2*k-1)*(-Ic.*CIsp + Is.*CIcp);
%  	denom = denom + coeff(2*k-1)*(-Ic.*CIs + Is.*CIc );
    
    numer = numer + coeff(2*k-1)*(Ic.*CIsp - Is.*CIcp);
 	denom = denom + coeff(2*k-1)*(Ic.*CIs - Is.*CIc );
    
%     numer = numer + coeff(2*k-2)*(Ic.*CIcp+Is.*CIsp) + coeff(2*k-1)*(-Ic.*CIsp + Is.*CIcp);
%  	denom = denom + coeff(2*k-2)*(Ic.*CIc +Is.*CIs ) + coeff(2*k-1)*(-Ic.*CIs + Is.*CIc );
%     numer(189, 232)
%     denom(189, 232)
%     r = numer(189, 232)/denom(189, 232)
end
g = numer./denom;
