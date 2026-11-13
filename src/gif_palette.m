function [map, lut] = gif_palette(S)
% GIF_PALETTE  common 256-colour palette for a set of frames (S = n x 3 uint8 sample pixels) and the look-up table lut
% common 256-colour palette for all frames: jet(112) (strain colours) + 32 greys + the most frequent colours of
% the sample frames; lut maps every 5-bit-per-channel RGB key (r*1024+g*32+b) to the nearest palette index
S = double(S);
fixed = [round(255*jet(112)); round(repmat(linspace(0,255,32).',1,3))];
key = floor(S(:,1)/8)*1024 + floor(S(:,2)/8)*32 + floor(S(:,3)/8);
cnt = accumarray(key+1, 1, [32768 1]);
[~,ord] = sort(cnt,'descend');
ord = ord(cnt(ord) > 0);
cols = zeros(0,3);
for q = 1:numel(ord)
    k = ord(q)-1;  c = [floor(k/1024), mod(floor(k/32),32), mod(k,32)]*8 + 4;
    if min(sum((fixed - c).^2,2)) > 3*8^2 && (isempty(cols) || min(sum((cols - c).^2,2)) > 3*8^2)
        cols(end+1,:) = c; %#ok<AGROW>
    end
    if size(fixed,1)+size(cols,1) >= 256, break, end
end
pal = [fixed; cols];
map = pal/255;
[kb,kg,kr] = ndgrid(0:31,0:31,0:31);                   % kb varies fastest
all5 = [kr(:) kg(:) kb(:)]*8 + 4;                       % row index = kr*1024+kg*32+kb+1
lut = zeros(32768,1,'uint8');
for a = 1:4096:32768
    b = min(a+4095,32768);
    D = zeros(b-a+1, size(pal,1));
    for c3 = 1:3
        D = D + (all5(a:b,c3) - pal(:,c3).').^2;
    end
    [~,im] = min(D,[],2);
    lut(a:b) = uint8(im-1);
end
end
