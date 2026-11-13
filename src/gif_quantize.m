function F = gif_quantize(A, lut)
% GIF_QUANTIZE  map a uint8 RGB frame to palette indices with the look-up table of gif_palette
key = bitshift(uint32(A(:,:,1)),-3)*1024 + bitshift(uint32(A(:,:,2)),-3)*32 + bitshift(uint32(A(:,:,3)),-3);
F = reshape(lut(double(key)+1), size(A,1), size(A,2));
end
