function run_real_snapshot(csvPath, pollutant, dateText)
% Conditionally run ONLY after verified ingestion: no fallback to synthetic.
% Sample usage: run_real_snapshot('data/processed/sima_utm.csv','PM2.5','2026-01-15')
% A single daily snapshot must be aggregated before using this routine.
if nargin<1,csvPath='data/processed/sima_utm.csv';end
if nargin<2,pollutant='PM2.5';end
if nargin<3,dateText='2026-01-15';end
if ~isfile(csvPath),error('NEXO:NoRealData','No authentic projected data at %s',csvPath);end
here=fileparts(mfilename('fullpath'));addpath(genpath(here));
t=readtable(csvPath,'TextType','string');
req={'station_id','timestamp_source','variable','value','unit','x_m','y_m'};
assert(all(ismember(req,t.Properties.VariableNames)), 'Input must use canonical metric CSV schema');
% Filtering by source calendar DATE avoids assuming the original zone is UTC.
match=contains(upper(t.variable),upper(string(pollutant))) & startsWith(t.timestamp_source,string(dateText));
t=t(match,:);
if isempty(t),error('NEXO:EmptySnapshot','No observations for requested pollutant and source date');end
assert(numel(unique(t.unit))==1,'Mixed units: reject without harmonization');
% Multiple hourly values per station within day are collapsed as explicit daily means.
[groups,stationLabels]=findgroups(t.station_id);
xs=splitapply(@mean,t.x_m,groups);ys=splitapply(@mean,t.y_m,groups);
y=splitapply(@(v)mean(v,'omitnan'),t.value,groups);
xy=[xs,ys];valid=all(isfinite(xy),2)&isfinite(y);xy=xy(valid,:);y=y(valid);
if numel(y)<6,error('NEXO:InsufficientStations','Fewer than six stations with valid observations');end
if size(unique(xy,'rows'),1)<6,error('NEXO:RepeatedCoords','Too few distinct station coordinates');end
[X,Y]=meshgrid(linspace(min(xs),max(xs),100),linspace(min(ys),max(ys),90));
q=[X(:),Y(:)];
% Width 15 km now expressed as 15000 m; do not mix m/km.
[z,~]=reconstruct(xy,y,q,'RBF',0.05,15000);Z=reshape(z,size(X));
h=convhull(xy(:,1),xy(:,2));inside=inpolygon(X,Y,xy(h,1),xy(h,2));Z(~inside)=NaN;
figure('Name','NEXO REAL DATA -- source day aggregate (verify data first)');
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],Z);set(gca,'YDir','normal');hold on;scatter(xy(:,1),xy(:,2),35,y,'filled','MarkerEdgeColor','w');hold off;axis equal tight;colorbar;
xlabel('UTM easting (m, EPSG:32614)');ylabel('UTM northing (m, EPSG:32614)');
title(sprintf('%s daily means; source calendar date %s -- observations pending QC',pollutant,dateText));
end
