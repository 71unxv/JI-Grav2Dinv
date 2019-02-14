function [g]=JIforwardGrav_blakely(model,dx,coor_m_x,coor_obs_x,dz,coor_m_z,coor_obs_z);

if length(dx)==1
    dx=ones(length(dx))*dx;
end
if length(dz)==1
    dz=ones(length(dz))*dz;
end

j=0;
for i=1:length(coor_x_obs);
kk=1;
progress_forward_grav=j*100/length(coor_x_obs)
    for j=1:length(coor_x_m);
        for k=1:length(coor_z_m);

            xk=coor_x_m(j);
            xobs=coor_x_obs(i);
            zk=coor_x_m(k);
            zobs=coor_z_obs(i);
            dxx=dx(j);
            dzz=dz(k);
            r=sqrt(((xk-xobs)^2)+((zk-zobs)^2));
g.kernel(i,kk)=2*yy*((zk-zobs)*dxx*dzz)/(r^2);
kk=kk+1;
        end
    end
    j=j+1;
end
g.rho=reshape(model,[(size(model,1))*(size(model,2)),1);
g.obs=g.kernel*g.obs;

end