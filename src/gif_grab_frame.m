function A = gif_grab_frame(fig)
% screen capture of the figure as uint8 RGB: getframe, else print to a temporary PNG
try
    F = getframe(fig); A = F.cdata;
    if isempty(A) || ndims(A) ~= 3, error('empty'); end
catch
    fn = [tempname() '.png'];
    print(fig, fn, '-dpng', '-r72');
    A = imread(fn); delete(fn);
end
A = uint8(A(:,:,1:3));
end
