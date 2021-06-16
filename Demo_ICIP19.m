clc; close all; clear all
% image = imread('4.bmp');
image = imread('5.png');

tic;
HSV = rgb2hsv( double(image) );   % RGB space to HSV  space
S = HSV(:,:,3);       % V layer
% [m,n]=size(S);
t0 = toc
f = S;

sigmas = 5;
w  = round(6*sigmas); if (mod(w,2) == 0); w  = w+1; end
T0  = 255; % maxFilter(f, w);

%% Input
sigmar = 30; eps = 1e-3;
T_set = [100:1:300];
samples = [-255:255]';
Err = 10* ones(1, length(T_set));
K = 1;
minErr = 10;
ii = sqrt(-1);

Kmax = 22;

while minErr > eps && K < Kmax
    
    Err = 10* ones(1, length(T_set));
    for j = 1:length(T_set)        
 
        w0 = 2*pi/(2*T_set(j) + 1);
        b =  exp(-0.5*samples.^2/sigmar^2).*sigmf(samples, [0.20, 0]); 
        A = ones(size(samples))./1;      

        for i = 2: K
            A = [A cos((i-1)*w0*samples) sin((i-1)*w0*samples)];
        end
        c = pinv(A)*b;
        recon = A*c;
        Err(j) = norm(recon - b, 2);
    end

    [minErr, Idx] = min(Err);
    K = K +1;
end

K = K-1;
T_star  = T_set(Idx)
N = 2*T_star +1;
w0 = 2*pi/N;
tic;
A = ones(size(samples));

for i = 2: K
    A = [A cos((i-1)*w0*samples) sin((i-1)*w0*samples)];
end



c = pinv(A)*b; 
t1 = toc
recon = A*c;

% plot(samples, recon, 'b', 'LineWidth', 2)
% hold on; grid on
% plot(samples, b, 'r', 'LineWidth', 2)
% legend('Reconstructed', 'Actual')
% title(['\sigma_r = ', num2str(sigmar)])
% err = b - recon;  
% max(abs(err))
% figure; stem(err)
% 


coeffnew = c';

tic
omega=(2*pi)/(2*T_set(Idx)+1);

g_opt = compress_Bright(f, coeffnew, sigmas, K, omega);
t2 = toc
% psnr(g, g_opt, 255)
L = g_opt;
L(L<0)=0;
L(L>255)=255;
% PSNR1 = psnr(g, L, 255)

% figure; imshow(uint8(g_opt))

tic;
L = g_opt;
R = S./L;

% Gamma correction

gamma = 2.2;
L_gamma = 255 * ( (L/255).^(1/gamma) );
enhanced_V = R .* L_gamma;
HSV(:,:,3) = enhanced_V;

f_our = hsv2rgb(HSV);  %  HSV space to RGB space
t5 = toc


 figure;
 imshow(uint8(image)),title('observed image');
 figure;
 imshow(uint8(f_our)),title('Proposed');


