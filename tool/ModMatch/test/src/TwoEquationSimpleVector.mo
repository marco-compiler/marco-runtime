package SimpleVector
	Real[10] vector(start = 10.0);
	equation
	for i in 1:5 loop
		vector[i] = 10;
	end for;
	for i in 1:5 loop
		vector[i+5] = 10;
	end for;
end SimpleVector;