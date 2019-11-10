% measure the saliency with spectral spatial distribution
function ssd = speSpaDis(img)
[path, ~] = fileparts(mfilename('fullpath'));
[nrow, ncol, nband] = size(img);
C = 4; % the number of endmember 
scale = 1/4; % reduce the image size to 1/4
im=imresize(img,scale); 
[nr,nc,nb] = size(im);
% VCA
vimg = reshape(im, [nr*nc, nb]);
M = vimg';
[ U, indicies, snrEstimate] = hyperVca(M, C);
flag = 1; % in case hyperFcls creats error. 
while(flag == 1)
try 
    [ AB ] = hyperFcls( M, U );
    flag = 0;
catch err
    [ U, indicies, snrEstimate] = hyperVca(M, C);  
    flag = 1;
end 
end
% L1/2  nmf
o.lambda = 0.5;
o.sources = C;
o.derta = 5;
o.S = AB;
nnschalfprior3(M, U, o, './testLhalf' );
load('./testLhalf');
U = A_est;
AB = S;
% f = figure; plot(U); 
% print(f, '-djpeg', 'Endmembers');
% title('VCA Recovered Endmembers'); grid on;
pcx = reshape(AB', [nr, nc, C]);
for i = 1:C
    imwrite(pcx(:,:,i), ['Abundance', num2str(i), '.jpg']);
end
mhc = 0;
for k=1:C
    Z(k) = sum(sum(pcx(:,:,k)));
end
m = nr; n = nc;
for k=1:C
    ind = find(ones(m,n)>0);
    [X, Y] = ind2sub([m,n], ind);
    pc = pcx(:,:,k);
    mh(k) = 1/Z(k)*sum(pc(:).*Y);
    vh(k) = 1/Z(k)*sum(pc(:).*(Y-mh(k)).*(Y-mh(k)));
    mv(k) = 1/Z(k)*sum(pc(:).*X);
    vv(k) = 1/Z(k)*sum(pc(:).*(X-mv(k)).*(X-mv(k)));
    D(k) = sum(pc(:).*sqrt(((X-m/2).*(X-m/2)+(Y-n/2).*(Y-n/2))));
end
V = vh+vv;
V = (V-min(V))/(max(V)-min(V));
D = (D-min(D))/(max(D)-min(D));

for i=1:m
    for j=1:n
        pc = pcx(i,j,:);
%         ssd(i,j) = sum(pc(:)'.*(1-V));%.*(1-D));
        ssd(i,j) = sum(pc(:)'.*(1-V).*(1-D));
    end
end
ssd = ssd/max(max(ssd));


