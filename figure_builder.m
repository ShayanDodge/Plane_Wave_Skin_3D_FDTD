clc
clear all
close all
uiopen('E:\University\project\Skin Modeling\1-derivative gaussian\dispersive\3-three_layer_dispersive_skin_fdtd_cpml_3D_gaussian\three_layer_dispersive_skin_fdtd_cpml_gaussian_3D_2D.fig',1)
for i=0:71
    saveas(gcf,['three_layer_dispersive_skin_fdtd_cpml_gaussian_3D_2D_' num2str(72-i) '.tiff'])
    close(gcf)
end