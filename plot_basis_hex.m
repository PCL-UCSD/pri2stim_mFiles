% plot_basis_hex
%
% assumes a square basis fcn

function [fh,ax] = plot_basis_hex(centerX,centerY,basis)

nf = max(max(abs(centerX)),max(abs(centerY)));
mv = min(min(centerX),min(centerY));
centerX = (0.75*centerX-mv)/(2*nf);
centerY = (0.75*centerY-mv)/(2*nf);

ax_size = [0.09 0.09];

centerX = centerX-ax_size(1)/2;centerY = centerY-ax_size(2)/2;

fh = figure;postmp = get(fh,'Position');set(fh,'Position',[postmp(1:3) postmp(3)]);%fh.Position(4) = fh.Position(3);

for ii = 1:length(centerX)
    
    ax(ii) = axes('Position',[centerX(ii) centerY(ii) ax_size]);
    %[centerX(ii) centerY(ii) ax_size]
    imagesc(reshape(basis(:,ii),sqrt(size(basis,1)),sqrt(size(basis,1))));
    axis square; axis xy; axis off;
    
    
end


return