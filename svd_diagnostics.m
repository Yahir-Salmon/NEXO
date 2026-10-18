function report=svd_diagnostics(A)
% Matrix diagnostics. QR-compatible MATLAB standard.
s=svd(A,'econ');report.s=s;report.rank=rank(A);
report.condition_number=cond(A);
report.singular_spectrum_ratio=s/max(s);
end
