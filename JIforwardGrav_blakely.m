function [g]=JIforwardGrav_blakely(model,dx,coor_x_m,coor_x_obs,dz,coor_z_m,coor_z_obs);
g.const=6.67e-11; %G;
if length(dx)==1;
    dx=ones(length(coor_x_m),1)*dx;
    y=1;
end
if length(dz)==1;
    dz=ones(length(coor_z_m),1)*dz;
end

jj=0;
for i=1:length(coor_x_obs);
kk=1;
progress_forward_grav=jj*100/length(coor_x_obs)
    for j=1:length(coor_x_m);
        for k=1:length(coor_z_m);
            xk=coor_x_m(j);
            xobs=coor_x_obs(i);
            zk=coor_z_m(k);
            zobs=coor_z_obs(i);
            dxx=dx(j);
            dzz=dz(k);
            r=sqrt(((xk-xobs)^2)+((zk-zobs)^2));
            g.kernel(i,kk)=2*g.const*((zk-zobs)*dxx*dzz)/(r^2);
            kk=kk+1;
        end
    end
    jj=jj+1;
end
g.rho=reshape(model,[(size(model,1))*(size(model,2)),1]);
g.obs=g.kernel*g.rho;
g.obs=g.obs*1e5;
end