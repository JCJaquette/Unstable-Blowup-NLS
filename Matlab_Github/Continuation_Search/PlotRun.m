x_grid= linspace(0,1,512);
x_len = length(x_grid);
t_len = 12568; %12568

Y=zeros(t_len ,x_len );
tic
parfor i =1:t_len 
    Y(i,:)=abs(u{i}(x_grid));
end
surf(Y)
colormap jet;
xlim([1,x_len ])
ylim([0,length(Y)])
view(2)
toc
colorbar
shading flat