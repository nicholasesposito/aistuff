subroutine eval_usplined(n,xs,p,q, x,y,dydx)
[eval_uspline]



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx

integer :: ia,ib
real(dp)::
xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm

if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; dydx=q(1); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; dydx=q(n); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
width
   if(xs(ib)>=x)exit 
found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2    
xr=x-xs(ia)-xh           
pm=(p(ia)+p(ib))*o2      
qm=(p(ib)-p(ia))/(xh*2) 
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
end subroutine eval_usplined


subroutine eval_usplinedd(n,xs,p,q, x,y,dydx,ddydxx)
[eval_uspline]



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx

integer :: ia,ib
real(dp)::
xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm

if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; dydx=q(1); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; dydx=q(n); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
width
   if(xs(ib)>=x)exit 
found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2    
xr=x-xs(ia)-xh           
pm=(p(ia)+p(ib))*o2      
qm=(p(ib)-p(ia))/(xh*2)  
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
ddydxx=qdh +qxh*xh*sh
end subroutine eval_usplinedd


subroutine eval_usplineddd(n,xs,p,q, x,y,dydx,ddydxx,dddydxxx)
[eval_uspline]



integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q
real(dp),             intent(in ):: x
real(dp),             intent(out):: y,dydx,ddydxx,dddydxxx

integer :: ia,ib
real(dp)::
xr,xh,pm,qm,qah,qbh,qxh,qdh,shh,chh,sh,ch,xcmsh,shm,chm,shhm,chhm

if(x<=xs(1))then; xr=x-xs(1); y=p(1)+q(1)*xr; dydx=q(1); return; endif
if(x>=xs(n))then; xr=x-xs(n); y=p(n)+q(n)*xr; dydx=q(n); return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
width
   if(xs(ib)>=x)exit 
found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2    
xr=x-xs(ia)-xh           
pm=(p(ia)+p(ib))*o2      
qm=(p(ib)-p(ia))/(xh*2)  
qah=q(ia)*o2;  qbh=q(ib)*o2
qxh=qah+qbh-qm 
qdh=qbh-qah    
shh=xh;        chh=u1
sh =xr;        ch =u1
shm =xr**3/6;  chm =xr**2*o2
shhm=xh**3/6;  chhm=xh**2*o2
xcmsh=xh**3/3
qdh=qdh/shh; qxh=qxh/xcmsh 
y=pm+xr*qm +qdh*(chm-chhm)  +  qxh*(xh*shm-xr*shhm)
dydx=qm+qdh*sh +qxh*(xh*chm-shhm)
ddydxx=qdh +qxh*xh*sh
dddydxxx=qxh*xh
end subroutine eval_usplineddd


subroutine eval_iuspline(n,xs, p,q,m,  x,y)
[eval_iuspline]




integer,              intent(in ):: n
real(dp),dimension(n),intent(in ):: xs,p,q,m
real(dp),             intent(in ):: x
real(dp),             intent(out):: y

real(dp),parameter:: u3o2=3*o2
real(dp):: a,b,c,d,t2,t3,t4,xh,xr,pa,pd,qa,qd
integer :: ia,ib

if(x<=xs(1))then; xr=x-xs(1); y=p(1)*xr+q(1)*xr**2/2; return; endif
if(x>=xs(n))then; xr=x-xs(n); y=m(n)+p(n)*xr+q(n)*xr**2/2; return; endif
do ib=2,n
   if(xs(ib)<=xs(ib-1))cycle 
width
   if(xs(ib)>=x)exit 
found
enddo
ia=ib-1
xh=(xs(ib)-xs(ia))*o2   
xr=x-xs(ia)-xh          
t2=xh**2/2
t3=t2*xh/3
pa=(p(ib)+p(ia))*o2
pd=(p(ib)-p(ia))*o2/xh
qa=(q(ib)+q(ia))*o2
qd=(q(ib)-q(ia))*o2/xh


c=qd
a=pa-c*t2
d=(qa-pd)*u3o2/t2
b=qa-d*t2
t2=xr**2/2
t3=t2*xr/3
t4=t3*xr/4
y=m(ia)+a*xr+b*t2+c*t3+d*t4
end subroutine eval_iuspline


subroutine best_tslalom(nh,mh,doru,hgts,hs,halfgate,bigT, & 
[best_slalom]
   hgtp,hp,qbest,yabest,enbest,modebest,maxita,maxitb,maxit,maxrts,FF)





integer,                 intent(in   ):: nh,mh,doru
integer, dimension(nh),  intent(in   ):: hgts
real(dp),dimension(nh),  intent(in   ):: hs
real(dp),                intent(in   ):: halfgate,bigT
integer, dimension(mh*2),intent(  out):: hgtp
real(dp),dimension(mh*2),intent(  out):: hp
real(dp),dimension(mh*2),intent(  out):: qbest
real(dp),dimension(mh*2),intent(  out):: yabest
real(dp),                intent(  out):: enbest
integer,dimension(mh),   intent(  out):: modebest
integer,                 intent(inout):: maxita,maxitb,maxit,maxrts
logical,                 intent(  out):: FF

integer, dimension(2,mh)  :: hgtn
real(dp),dimension(mh*2)  :: q,ya
real(dp),dimension(2,2,mh):: hn
real(dp)                  :: en,tspan,hspan,enbase,hgbigT
integer, dimension(mh)    :: code,mode
integer, dimension(mh*2)  :: bend
integer                   :: i,k,m,route_count,ita,ittot
logical, dimension(mh*2)  :: off
logical                   :: flag,descending

m=mh*2
call set_gates(nh,mh,doru,hgts,hs, hgtn,hn,code,FF)


if    (hn(1,2,1)>hn(1,1,mh))then; descending=T 
elseif(hn(2,2,1)<hn(2,1,mh))then; descending=F 
else 
   descending=(doru==1)
endif
hgbigT=bigT/halfgate 
tspan=(hgtn(2,mh)-hgtn(1,1))/hgbigT
if(descending)then; hspan=(hn(1,1,1)-hn(1,2,mh))
else              ; hspan=(hn(2,2,mh)-hn(2,1,1))
endif
enbase=enbase_t(tspan,hspan) 
if(FF)then
   write(41,*) 'In best_tslalom; failure flag was raised in call to
set_gates'
   return
endif
call count_routes(mh,code,route_count,FF)
maxrts=max(maxrts,route_count)
if(FF)then
   write(41,*)&
        'In best_tslalom; failure flag was raised in call to
count_routes'
   return
endif
if(route_count>4)call list_routes(mh,code) 
when >4
enbest=hu
flag=T
do k=1,ihu
   call next_route(mh,code,mode,flag)
   if(flag)then; flag=F; exit; endif
   call set_posts(mh,mode,hgtn,hn,bend,hgtp,hp,off)
   call slalom_tspline(m,bend,hgtp,hp,off,hgbigT, &
        q,ya,en,ita,maxitb,ittot,FF); en=en/enbase
   maxita=max(maxita,ita)
   maxit =max(maxit,ittot)
   if(FF)then
      write(41,*) &
           'In best_tslalom; failure flag was raised in call to
slalom_tspline'
      return
   endif
   if(en<enbest)then
      modebest=mode
      enbest  =en
      qbest   =q
      yabest  =ya
   endif
enddo
end subroutine best_tslalom

subroutine best_uslalom(nh,mh,doru,hgts,hs,halfgate,  & 
     hgtp,hp,qbest,yabest,enbest,modebest,maxita,maxitb,maxit,maxrts,FF)




integer,                 intent(in   ):: nh,mh,doru
integer, dimension(nh),  intent(in   ):: hgts
real(dp),dimension(nh),  intent(in   ):: hs
real(dp),                intent(in   ):: halfgate
integer, dimension(mh*2),intent(  out):: hgtp
real(dp),dimension(mh*2),intent(  out):: hp
real(dp),dimension(mh*2),intent(  out):: qbest
real(dp),dimension(mh*2),intent(  out):: yabest
real(dp),                intent(  out):: enbest
integer,dimension(mh),   intent(  out):: modebest
integer,                 intent(inout):: maxita,maxitb,maxit,maxrts
logical,                 intent(  out):: FF

integer, dimension(2,mh)  :: hgtn
real(dp),dimension(mh*2)  :: q,ya
real(dp),dimension(2,2,mh):: hn
real(dp)                  :: en
integer, dimension(mh)    :: code,mode
integer, dimension(mh*2)  :: bend
integer                   :: i,k,m,route_count,ita,ittot
logical, dimension(mh*2)  :: off
logical                   :: flag

m=mh*2
call set_gates(nh,mh,doru,hgts,hs, hgtn,hn,code,FF)
if(FF)then
   write(41,*) 'In best_uslalom; failure flag was raised in call to set_gates'
   return
endif
call count_routes(mh,code,route_count,FF)
maxrts=max(maxrts,route_count)
if(FF)then
   write(41,*)&
        'In best_uslalom; failure flag was raised in call to count_routes'
   return
endif
if(route_count>4)call list_routes(mh,code)
enbest=hu
flag=T
do k=1,ihu
   call next_route(mh,code,mode,flag)
   if(flag)then; flag=F; exit; endif
   call set_posts(mh,mode,hgtn,hn,bend,hgtp,hp,off)
   call slalom_uspline(m,bend,hgtp,hp,off,halfgate, q,ya,en,ita,maxitb,ittot,FF)
   maxita=max(maxita,ita)
   maxit =max(maxit,ittot)
   if(FF)then
      write(41,*) &
           'In best_uslalom; failure flag was raised in call to slalom_uspline'
      return
   endif
   if(en<enbest)then
      modebest=mode
      enbest  =en
      qbest   =q
      yabest  =ya
   endif
enddo
end subroutine best_uslalom


subroutine count_gates(nh,hgts,mh)




integer,              intent(in ):: nh
integer,dimension(nh),intent(in ):: hgts
integer,              intent(out):: mh

integer:: hgtp
integer:: i

hgtp=hgts(1)-1 
mh=0
do i=1,nh
   if(hgts(i)<=hgtp)cycle

   mh=mh+1
   hgtp=hgts(i)
enddo
end subroutine count_gates
