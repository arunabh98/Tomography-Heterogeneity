function A = assembleA(thetas, order)
    
	A = zeros(order * length(thetas)  ,((order+1)*(order+2)/2) - 1);
	colstart = 1;
	for i = 1:order
		Ai = zeros(length(thetas), (i+1));
		for j = 0:i
			Ai(:,(j+1)) = nchoosek(i,j) * (sind(thetas)).^j .* (cosd(thetas)).^(i-j);
		end
		A( (((i-1)*length(thetas))+1) : (i*length(thetas)) , colstart:(colstart+size(Ai,2)-1) ) = Ai;
		colstart = colstart + size(Ai,2);
	end
    
end
    
