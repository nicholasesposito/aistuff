subroutine set_gates(nh,mh,doru,hgts,hs, hgtn,hn,code,FF)
integer,                   intent(in ):: nh,mh,doru
integer, dimension(nh),    intent(in ):: hgts
real(dp),dimension(nh),    intent(in ):: hs
integer, dimension(2,  mh),intent(out):: hgtn
real(dp),dimension(2,2,mh),intent(out):: hn
integer, dimension(    mh),intent(out):: code
logical,                   intent(out):: FF
real(dp):: hp
integer :: i,im,i2,i2m,imh,n,atti,attim,codeim,hgtp
FF=F
n=nh*2
hgtp=hgts(1)-1 
imh=0
do i=1,nh
   i2=i*2
   i2m=i2-1
   hp=hs(i)
   if(hgts(i)>hgtp)then
      imh=imh+1
      hgtp=hgts(i)
      hgtn(1,imh)=hgtp-1
      hgtn(2,imh)=hgtp+1
      hn(:,:,imh)=hp
   elseif(hgts(i)<hgtp) then
      FF=T
      write(41,*) 'In set_gates; data are not temporally monotonic'
      return
   else
      hn(1,1,imh)=max(hn(1,1,imh),hp)
      hn(1,2,imh)=min(hn(1,2,imh),hp)
      hn(2,2,imh)=max(hn(2,2,imh),hp)
      hn(2,1,imh)=min(hn(2,1,imh),hp)
   endif
enddo
if(imh/=mh)stop 'In set_gates; inconsistent gate tallies, imh and mh'
if(mh==1)then
   code(1)=4*doru
   return
endif
attim=0
code=0 
codeim=9 
do i=2,mh
   atti=0 
   im=i-1
   if(hgtn(1,i)<=hgtn(2,im))then
      if(hn(2,2,im)<=hn(1,2,i))then
         atti=2 
         code(i)=2
         if(attim==2.and.(codeim==0.or.codeim==2))code(im)=8
      elseif(hn(2,1,im)>=hn(1,1,i))then
         atti=1 
         code(i)=3
         if(attim==1.and.(codeim==0.or.codeim==3))code(im)=4
      else
         code(i)=5
         if(hn(2,1,im)<=hn(1,2,i))then; hn(1,2,i) =hn(2,1,im)
                                  else; hn(2,1,im)=hn(1,2,i)
         endif
         if(hn(2,2,im)<=hn(1,1,i))then; hn(2,2,im)=hn(1,1,i)
                                  else; hn(1,1,i) =hn(2,2,im)
         endif
      endif
   else
      if(hn(2,2,im)<=hn(1,2,i))then
         atti=2 
         if(attim==2.and.(codeim==0.or.codeim==2))code(im)=8
      elseif(hn(2,1,im)>=hn(1,1,i))then
         atti=1 
         if(attim==1.and.(codeim==0.or.codeim==3))code(im)=4
      endif
   endif
   attim=atti
   codeim=code(i)
enddo
end subroutine set_gates


subroutine set_posts(mh,mode,hgtn,hn, bend,hgtp,hp,off)
[set_posts]
integer,                   intent(in ):: mh
integer, dimension(    mh),intent(in ):: mode
integer, dimension(2,  mh),intent(in ):: hgtn
real(dp),dimension(2,2,mh),intent(in ):: hn
integer, dimension(mh*2),  intent(out):: bend,hgtp
real(dp),dimension(mh*2),  intent(out):: hp
logical, dimension(mh*2),  intent(out):: off
real(dp):: hprev
integer :: i,i2,i2m,i2mm,im,modei,hgtprev
off=F
do i=1,mh
   im=i-1
   modei=mode(i)
   i2=i*2; i2m=i2-1; i2mm=i2-2
   hgtp(i2m)=hgtn(1,i)
   hgtp(i2 )=hgtn(2,i)
   hp(i2m)=hn(1,modei,i)
   hp(i2 )=hn(2,modei,i)
   if(i>1)then
      if(hgtprev==hgtp(i2m))then
         if(hprev==hp(i2m))off(i2m)=T
         if(mode(im)==2.and.modei==1)then
            if(hprev<=hp(i2m))then
               off(i2mm)=T
            else
               off(i2m)=T
            endif
         elseif(mode(im)==1.and.modei==2)then
            if(hprev<=hp(i2m))then
               off(i2m)=T
            else
               off(i2mm)=T
            endif
         endif
      endif
   endif
   bend(i2m)=modei*2-3 
   bend(i2 )=-bend(i2m)
   hgtprev=hgtp(i2)
   hprev  =hp(i2)
enddo
end subroutine set_posts


subroutine count_routes(n,code,count,FF)
integer,             intent(in ):: n
integer,dimension(n),intent(in ):: code
integer,             intent(out):: count
logical,             intent(out):: FF
integer,dimension(n):: mode
logical             :: flag
FF=F
flag=T
do count=0,ihu; call next_route(n,code,mode,flag); if(flag)return; enddo
FF=(count>ihu)
if(FF) write(41,*) 'In count_routes; number of routes exceeds allowance
= ',ihu
end subroutine count_routes


subroutine list_routes(n,code)
integer,             intent(in ):: n
integer,dimension(n),intent(in ):: code
integer,dimension(n):: mode
integer             :: i
logical             :: flag
write(41,'("List all route combinations of ",i4," allowed passage
modes")'),n
flag=T
do i=1,ihu
   call next_route(n,code,mode,flag)
   if(flag)then
      write(41,'(" In list_routes; List of routes complete")'); flag=F;
exit
   endif
   write(41,60)i,mode
enddo
if(i>ihu) write(41,'("This list is not necessarily complete")')
60 format(i5,3x,6(2x,5i2))
end subroutine list_routes


subroutine next_route(n,code,mode,flag)
integer,             intent(in   ):: n
integer,dimension(n),intent(in   ):: code
integer,dimension(n),intent(inout):: mode
logical,             intent(inout):: flag
integer,dimension(0:8,2):: options 
code
integer,dimension(0:2)  :: firstmode
integer                 :: i,im,j,modeim,modejm,option
data options/0,1,2,0,1,2,0,1,2, 0,0,0,1,1,1,2,2,2/
data firstmode/1,1,2/
modeim=1 
if(flag)then
   do i=1,n
      option=options(code(i),modeim)
      mode(i)=firstmode(option)
      modeim=mode(i)
   enddo
   flag=F
   return
endif
do i=n,1,-1
   im=i-1
   if(i>1)then
      modeim=mode(im)
   else
      modeim=1
   endif
   option=options(code(i),modeim)
   if(option>0.or.mode(i)==2)cycle
   mode(i)=2
   modejm=mode(i)
   do j=i+1,n
      option=options(code(j),modejm)
      mode(j)=firstmode(option)
      modejm=mode(j)
   enddo
   return
enddo
flag=T
end subroutine next_route


subroutine slalom_tspline(n,bend,hgxn,yn,off,bigX, &
     q,ya,en,ita,maxitb,ittot,FF)                               
integer,                intent(in   ):: n
integer, dimension(n),  intent(in   ):: bend,hgxn
real(dp),dimension(n),  intent(in   ):: yn
logical, dimension(n),  intent(in   ):: off
real(dp),               intent(in   ):: bigX
real(dp),dimension(n),  intent(  out):: q
real(dp),dimension(n),  intent(  out):: ya
real(dp),               intent(  out):: en
integer,                intent(  out):: ita,ittot
integer,                intent(inout):: maxitb
logical,                intent(  out):: FF
integer,parameter      :: nita=50,nitb=80
real(dp),dimension(n)  :: xs,jump,qt,yat
real(dp)               :: sj,sjmin,ena
integer                :: i,j,k,itb,hgxp
logical,dimension(n)   :: on
FF=F
xs=hgxn/bigX
hgxp=hgxn(1)-1
do i=1,n
   if(off(i))then; on(i)=F; cycle; endif
   on(i)=(hgxn(i)>hgxp); if(on(i))hgxp=hgxn(i)
enddo
ittot=1
call fit_gtspline(n,xs,yn,on,qt,jump,yat,en,FF)
ena=en
if(FF)then
   write(41,*) 'In slalom_tspline; failure flag raised in call to
fit_gtspline'
   write(41,*) 'at initialization of A loop'
   return
endif
do ita=1,nita
   q=qt   
   ya=yat 
   j=0
   k=0
   sjmin=0
   do i=1,n
      if(.not.on(i))cycle
      sj=-bend(i)*jump(i)
      if(sj<0)then
         j=i
         on(i)=F
      else
         k=k+1 
      endif
   enddo
   if(j==0)exit 
   if(k==0)on(j)=T 
   do itb=1,nitb
      call fit_gtspline(n,xs,yn,on,qt,jump,yat,en,FF)
      if(FF)then
         write(41,*)&
              'In slalom_tspline; failure flag raised in call to
fit_gtspline'
         write(41,*) 'at B loop, iterations ita,itb = ',ita,itb
         return
      endif
       ittot=ittot+1 
fit_tspline
      j=0
      sjmin=u1
      do i=1,n
         if(on(i).or.off(i))cycle
         sj=bend(i)*(yn(i)-yat(i))
         if(sj<0)then
            sj=(yn(i)-ya(i))/(yat(i)-ya(i))
            if(sj<sjmin)then
               j=i
               sjmin=sj
            endif
         endif
      enddo
      if(j==0)exit 
solution as A
      ya=ya+sjmin*(yat-ya)
      q=q+sjmin*(qt-q)
      on(j)=T
   enddo 
   maxitb=max(maxitb,itb)
   if(itb>nitb) then
      FF=T
      write(41,*) 'In slalom_tspline; exceeding the allocation of B
iterations'
      return
   end if
   q=qt
   ya=yat
   if(en>=ena)then
      write(41,*) 'In slalom_tspline; energy failed to decrease'
      exit 
   endif
   ena=en
enddo 
if(ita>nita)then
   FF=T
   write(41,*) 'In slalom_tspline; exceeding the allocation of A
iterations'
   return
endif
end subroutine slalom_tspline


