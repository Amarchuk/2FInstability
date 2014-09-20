;MAJOR AXIS################################################################################################
s=read_asc('resNGC3245pa355_stars_sn25.txt',1)
g=read_asc('resNGC3245pa355_gas.txt',1)
rsec=s[0,*]
vel=s[2,*]
evel=s[3,*]
sig=s[4,*]
esig=s[5,*]
rgsec=g[0,*]
gvhb=g[1,*]
egvhb=g[2,*]
gvoi=g[7,*]
egvoi=g[8,*]
set_plot,'PS'
device,filename='fig_n3245_major.ps',/col,xs=15,ys=15
!x.thick=3
!y.thick=3
!p.thick=3
!p.font=0
!p.multi=[0,1,2]
!p.charsize=1.15
xra=[-70,70]
xtint=20
vsys=1300
yr=300
  plotsym,0,0.5,/fill
  plot,rsec,vel,/nodata,yr=[vsys-yr,vsys+yr],xr=xra,xtickint=xtint,xst=1,yst=1,xtickformat='(A1)',$
       ymargin=[1.75,1],ytit='Velocity, km/s'
    oplot,rsec,vel,ps=8
    oploterror,rsec,vel,evel,psym=3,hatlen=80
  plotsym,4,0.7
    oplot,rgsec,gvhb,ps=8
    oploterror,rgsec,gvhb,egvhb,psym=8,hatlen=80
  plotsym,8,0.7
    oplot,rgsec,gvoi,ps=8
    oploterror,rgsec,gvoi,egvoi,psym=8,hatlen=80
    oplot,[0,0],vsys+[-500,500],linest=2,thick=1
    oplot,[-100,100],[vsys,vsys],linest=2,thick=1
  plotsym,0,0.8,/fill
  legend,['stars','H!9b!3','[OIII]'],psym=[8,5,6],/bottom,box=0

  plot,rsec,sig,/nodata,yr=[0,230],xr=xra,xtickint=xtint,xst=1,yst=1,ymargin=[3.5,-0.75],xtit='Radial distances, arcsec',$
      ytit='Dispersion, km/s'
  plotsym,0,0.5,/fill
    oplot,rsec,sig,ps=8
    oploterror,rsec,sig,esig,psym=3,hatlen=100
  plotsym,0.0,0.5
    oplot,[0,0],[0,300],linest=2,thick=1

!p.multi=0
device,/close
;MINOR AXIS################################################################################################
s=read_asc('resNGC3245pa267_stars_sn25.txt',1)
g=read_asc('resNGC3245pa267_gas.txt',1)
cor=20
rsec=s[0,*]
vel=s[2,*]+cor
evel=s[3,*]
sig=s[4,*]
esig=s[5,*]
rgsec=g[0,*]
gvhb=g[1,*]+cor+70
egvhb=g[2,*]
gvoi=g[7,*]+cor
egvoi=g[8,*]
set_plot,'PS'
device,filename='fig_n3245_minor.ps',/col,xs=15,ys=15
!x.thick=3
!y.thick=3
!p.thick=3
!p.font=0
!p.multi=[0,1,2]
xra=[-40,40]
xtint=10
vsys=1300
yr=300
  plotsym,0,0.5,/fill
  plot,rsec,vel,/nodata,yr=[vsys-yr,vsys+yr],xr=xra,xtickint=xtint,xst=1,yst=1,xtickformat='(A1)',$
       ymargin=[1.75,1],ytit='Velocity, km/s'
    oplot,rsec,vel,ps=8
    oploterror,rsec,vel,evel,psym=3,hatlen=80
  plotsym,4,0.7
    oplot,rgsec,gvhb,ps=8
    oploterror,rgsec,gvhb,egvhb,psym=8,hatlen=80
  plotsym,8,0.7
    oplot,rgsec,gvoi,ps=8
    oploterror,rgsec,gvoi,egvoi,psym=8,hatlen=80
    oplot,[0,0],vsys+[-500,500],linest=2,thick=1
    oplot,[-100,100],[vsys,vsys],linest=2,thick=1
  plotsym,0,0.8,/fill
  legend,['stars','H!9b!3','[OIII]'],psym=[8,5,6],/bottom,box=0

  plot,rsec,sig,/nodata,yr=[0,230],xr=xra,xtickint=xtint,xst=1,yst=1,ymargin=[3.5,-0.75],xtit='Radial distances, arcsec',$
      ytit='Dispersion, km/s'
  plotsym,0,0.5,/fill
    oplot,rsec,sig,ps=8
    oploterror,rsec,sig,esig,psym=3,hatlen=100
  plotsym,0.0,0.5
    oplot,[0,0],[0,300],linest=2,thick=1
!p.multi=0
device,/close
set_plot,'X'
end
