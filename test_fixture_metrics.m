% MATLAB equivalent fixture parity, execute from project root.
a=readtable('data/synthetic/metric_fixture.csv');
ref=readtable('data/synthetic/metric_fixture_expected.csv');
diff=a.observed-a.predicted;
assert(abs(sqrt(mean(diff.^2))-ref.rmse)<1e-10);
assert(abs(mean(abs(diff))-ref.mae)<1e-10);
fprintf('Fixture metric parity: MATLAB matches Python reference.\n');
