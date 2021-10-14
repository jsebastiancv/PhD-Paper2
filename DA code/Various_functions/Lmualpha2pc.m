function pc = Lmualpha2pc(L, mu, alpha)
%LMUALPHA2PC convert L, mu and alpha to pc
	mc2 = 0.511;
	pc = sqrt(2.0 .* mu .* mc2 .* B0(L)) ./ sin(alpha);
return
