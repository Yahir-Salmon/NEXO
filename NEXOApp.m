function NEXOApp()
% NEXO scientific UI. Requires MATLAB R2020b+ for uifigure/UI controls.
% No commercial toolboxes. ALL current rendered sample points are SYNTHETIC.
% Start: addpath(genpath('matlab')); NEXOApp
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
addpath(genpath(fileparts(here)));
fig=uifigure('Name','NEXO | Laboratorio de reconstrucción ambiental - SINTÉTICO',...
    'Position',[100 80 1350 780],'Color',[0.07 0.10 0.15]);
main=uigridlayout(fig,[1,2]);main.ColumnWidth={265,'1x'};main.RowHeight={'1x'};main.Padding=[12 12 12 12];
ctrl=uipanel(main,'Title','NEXO · Parámetros reproducibles');ctrl.Layout.Column=1;
g=uigridlayout(ctrl,[17 1]);g.RowHeight={30,40,23,38,23,38,23,38,23,38,23,38,35,35,35,'1x',30};
uiwrap=@(s) uilabel(g,'Text',s,'FontWeight','bold');
uiwrap('Modo: laboratorio sintético');
status=uilabel(g,'Text','SIN DATOS REALES. Semilla fija.','FontColor',[.60 .26 .12]);
uiwrap('Número de estaciones (mínimo 8)');
ns=uispinner(g,'Limits',[8,60],'RoundFractionalValues','on','Step',1,'Value',24);
uiwrap('Ruido de medición (µg/m³)');
noise=uislider(g,'Limits',[0,8],'Value',1,'MajorTicks',[0 2 4 6 8]);
uiwrap('Regularización Tikhonov, λ');
lam=uidropdown(g,'Items',{'0.0001','0.001','0.01','0.05','0.1','0.3','1','3','10'},'Value','0.05');
uiwrap('Reconstrucción');
method=uidropdown(g,'Items',{'RBF','IDW'},'Value','RBF');
uiwrap('Velocidad u (m/s)');
wind=uislider(g,'Limits',[-4 4],'Value',1.5,'MajorTicks',-4:2:4);
uiwrap('Difusión D (m²/s)');
diffus=uislider(g,'Limits',[0 60],'Value',10,'MajorTicks',0:20:60);
runbtn=uibutton(g,'Text','Ejecutar experimento','ButtonPushedFcn',@refresh,'FontWeight','bold');
uibutton(g,'Text','Exportar figura PNG','ButtonPushedFcn',@export_current);
uibutton(g,'Text','Guardar escenario JSON','ButtonPushedFcn',@save_scenario);
uibutton(g,'Text','Cargar escenario JSON','ButtonPushedFcn',@load_scenario);
seedlbl=uilabel(g,'Text','Semilla reproducible: 20261018');
right=uitabgroup(main);right.Layout.Column=2;
tabHome=uitab(right,'Title','Inicio / datos');
tabMap=uitab(right,'Title','Reconstrucción');
tabSvd=uitab(right,'Title','Álgebra SVD');
tabGrad=uitab(right,'Title','Cálculo vectorial');
tabPde=uitab(right,'Title','Transporte-difusión');
tabCV=uitab(right,'Title','Validación');
tabReport=uitab(right,'Title','Informe / fuentes');
axMap=uiaxes(tabMap,'Position',[50 75 925 585]);axSvd=uiaxes(tabSvd,'Position',[50 75 925 585]);
axGrad=uiaxes(tabGrad,'Position',[50 75 925 585]);axPde=uiaxes(tabPde,'Position',[50 75 925 585]);
axCV=uiaxes(tabCV,'Position',[50 75 925 585]);
intro=uitextarea(tabHome,'Editable','off','Position',[30 65 900 600],...
    'Value',{'NEXO | Laboratorio de reconstrucción matemática y modelación ambiental urbana';'';...
    'VERSION 0.1: demostrador científico sobre función CONOCIDA Y SINTÉTICA.';...
    'Las ubicaciones de los puntos son simuladas en coordenadas locales de kilómetros.';...
    'No son estaciones del SIMA ni se infieren concentraciones reales en Monterrey.';'';...
    'Reproducible: 24 puntos aleatorios, semilla 20261018, ruido gaussiano.';...
    'Métodos: IDW y RBF gaussiana con Tikhonov por SVD.';...
    'Pruebas: dejando una estación fuera; error contra función sintética conocida.';...
    'EDP: transporte por upwind conservativo, difusión centrada y fronteras periódicas.';...
    'Ajusta controles a la izquierda y pulsa Ejecutar; cambia de pestaña.';'';...
    'CONEXIÓN A DATOS REALES: scripts de descarga independientes en src/acquisition/.';...
    'Descarga y validación deben completarse ANTES de afirmar resultados sobre Monterrey.'});
uitextarea(tabReport,'Editable','off','Position',[25 75 940 570],...
    'Value',{'NEXO | DOCUMENTACIÓN Y FUENTES';'';...
    'INECC: datos.gob.mx - concentraciones validadas en Nuevo León (CSV, 2022).';...
    'Observatorio del Aire: observatoriodelaire.com/datos - PM2.5/PM10 2025-2026.';...
    'NASA POWER: power.larc.nasa.gov - meteorología de reanálisis/modelo.';...
    'OpenAQ v3 requiere credencial API (X-API-Key).';'';...
    'Esta app no muestra datos medidos. No extrapolar al área metropolitana.';...
    'El informe editable se encuentra en paper/main.tex.';...
    'Los resultados ejecutados en Python, que no reemplazan MATLAB, están en results/.'});
right.SelectionChangedFcn=@refresh;
refresh();

    function refresh(~,~)
        dlg=[];
        try
            dlg=uiprogressdlg(fig,'Title','NEXO','Message','Reconstrucción y validación...',...
                'Indeterminate','on');
            n=round(ns.Value);sd=noise.Value;reg=str2double(lam.Value);
            [xy,y]=synthetic_stations(n,20261018,sd,false);
            [X,Y]=meshgrid(linspace(1,59,92),linspace(1,49,78));q=[X(:),Y(:)];
            [est,~]=reconstruct(xy,y,q,method.Value,reg,15);
            Z=reshape(est,size(X));
            % No convex hull rendering outside the sampled support.
            try
                hull=convhull(xy(:,1),xy(:,2));
                in=inpolygon(X,Y,xy(hull,1),xy(hull,2));
            catch, in=true(size(X));end
            Z(~in)=NaN;
            cla(axMap);imagesc(axMap,[X(1,1),X(1,end)],[Y(1,1),Y(end,1)],Z);
            set(axMap,'YDir','normal');hold(axMap,'on');scatter(axMap,xy(:,1),xy(:,2),40,y,'filled','MarkerEdgeColor','white');hold(axMap,'off');
            colorbar(axMap);axis(axMap,'equal');axis(axMap,'tight');xlabel(axMap,'x local (km)');ylabel(axMap,'y local (km)');
            title(axMap,['SINTÉTICO · ' method.Value ' · fuera de la envolvente: sin estimación']);
            A=rbf_design(xy,xy,15);sp=svd(A);
            cla(axSvd);semilogy(axSvd,1:numel(sp),sp,'o-','LineWidth',1.4);grid(axSvd,'on');
            xlabel(axSvd,'índice singular');ylabel(axSvd,'σᵢ (adimensional)');
            title(axSvd,sprintf('SVD del diseño RBF | cond(A)=%.3g | rango=%d',cond(A),rank(A)));
            % Gradients expressed per km (NOT per m).
            [Gx,Gy,mag]=field_gradient(fillmissing(Z,'nearest'),X(1,2)-X(1,1),Y(2,1)-Y(1,1));
            mag(~in)=NaN;cla(axGrad);imagesc(axGrad,[X(1,1),X(1,end)],[Y(1,1),Y(end,1)],mag);set(axGrad,'YDir','normal');
            hold(axGrad,'on');skip=7;quiver(axGrad,X(1:skip:end,1:skip:end),Y(1:skip:end,1:skip:end),...
                Gx(1:skip:end,1:skip:end),Gy(1:skip:end,1:skip:end),'k');hold(axGrad,'off');
            axis(axGrad,'equal');axis(axGrad,'tight');colorbar(axGrad);
            title(axGrad,'SINTÉTICO · |grad C| [µg/m³/km] · derivado del interpolante');
            xlabel(axGrad,'x local (km)');ylabel(axGrad,'y local (km)');
            pde=run_periodic_experiment(64,wind.Value,diffus.Value,500);
            cla(axPde);imagesc(axPde,pde.C);set(axPde,'YDir','normal');colorbar(axPde);
            title(axPde,sprintf('EDP SINTÉTICA | masa relativa=%.2g | error RMSE=%.3f µg/m³',...
                pde.relative_mass_error,pde.rmse));xlabel(axPde,'índice de celda x');ylabel(axPde,'índice de celda y');
            cv1=loocv_stations(xy,y,'IDW',reg,15);
            cv2=loocv_stations(xy,y,'RBF',reg,15);
            cla(axCV);bar(axCV,[cv1.rmse,cv2.rmse]);
            set(axCV,'XTick',1:2,'XTickLabel',{'IDW','RBF + Tikhonov'});ylabel(axCV,'RMSE LOOCV (µg/m³)');
            title(axCV,'SINTÉTICO · predicción de una estación excluida');grid(axCV,'on');
            status.Text=sprintf('n=%d | λ=%.4g | RMSE LOOCV=%.2f',n,reg,strcmp(method.Value,'IDW')*cv1.rmse+strcmp(method.Value,'RBF')*cv2.rmse);
            if ~isempty(dlg),close(dlg);end
        catch ME
            if ~isempty(dlg),close(dlg);end
            uialert(fig,ME.message,'Error de cálculo');
        end
    end
    function export_current(~,~)
        t=right.SelectedTab;ax=axMap;
        if t==tabSvd,ax=axSvd;elseif t==tabGrad,ax=axGrad;elseif t==tabPde,ax=axPde;elseif t==tabCV,ax=axCV;end
        [f,p]=uiputfile('*.png','Exportar gráfica sintética');
        if isequal(f,0),return;end
        exportgraphics(ax,fullfile(p,f),'Resolution',200);
    end
    function save_scenario(~,~)
        [f,p]=uiputfile('*.json','Guardar configuración de experimento');
        if isequal(f,0),return;end
        scenario=struct('kind','SYNTHETIC','seed',20261018,'n',ns.Value,...
            'noise_std',noise.Value,'lambda',str2double(lam.Value),...
            'method',method.Value,'wind_u',wind.Value,'diffusivity',diffus.Value);
        fid=fopen(fullfile(p,f),'w');fprintf(fid,'%s',jsonencode(scenario));fclose(fid);
    end
    function load_scenario(~,~)
        [f,p]=uigetfile('*.json','Cargar escenario');if isequal(f,0),return;end
        try
            scenario=jsondecode(fileread(fullfile(p,f)));
            assert(strcmp(scenario.kind,'SYNTHETIC'),'Unsupported scenario kind');
            ns.Value=scenario.n;noise.Value=scenario.noise_std;
            lam.Value=num2str(scenario.lambda);method.Value=scenario.method;
            wind.Value=scenario.wind_u;diffus.Value=scenario.diffusivity;
            refresh();
        catch ME,uialert(fig,ME.message,'Escenario inválido');end
    end
end
