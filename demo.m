scene = importdata('scene18.mat');
img = imread('scene18_RGB.jpg');
[HSI_ssd,Sm_HSI_EAS] = Saliency_EAS(scene);
imgS = imresize(img, size(Sm_HSI_EAS));
figure, 
imshowpair(imgS,Sm_HSI_EAS, 'montage');