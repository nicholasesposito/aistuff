subroutine slalom_uspline(n,bend,hgxn,yn,off,halfgate,&
     q, ya,en,ita,maxitb,ittot,FF)       
integer,                intent(in   ):: n
integer, dimension(n),  intent(in   ):: bend,hgxn
real(dp),dimension(n),  intent(in   ):: yn
logical, dimension(n),  intent(in   ):: off
real(dp),               intent(in   ):: halfgate
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
xs=hgxn*halfgate
hgxp=hgxn(1)-1
do i=1,n
   if(off(i))then
      on(i)=F
      cycle
   endif
   on(i)=(hgxn(i)>hgxp)
   if(on(i))hgxp=hgxn(i)
enddo
ittot=1
call fit_guspline(n,xs,yn,on,qt,jump,yat,en,FF)
ena=en
if(FF)then
   write(41,*) 'In slalom_uspline; failure flag raised in call to
fit_guspline'
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
      call fit_guspline(n,xs,yn,on,qt,jump,yat,en,FF)
      if(FF)then
         write(41,*)&
              'In slalom_uspline; failure flag raised in call to
fit_guspline'
         write(41,*) 'at B loop, iterations ita,itb = ',ita,itb
         return
      endif
      ittot=ittot+1 
fit_uspline
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
      write(41,*) 'In slalom_uspline; exceeding the allocation of B
iterations'
      return
   end if
   q=qt
   ya=yat
   if(en>=ena)then
      write(41,*) 'In slalom_uspline; energy failed to decrease'
      exit
   endif
   ena=en
enddo
if(ita>nita)then
   FF=T
   write(41,*) 'In slalom_uspline; exceeding the allocation of A
iterations'
   return
endif
end subroutine slalom_uspline


subroutine convertd(n,halfgate,tdata,hdata,phof,&
     doru,idx,hgts,hs,descending,FF)
integer,              intent(in ):: n
real(dp),             intent(in ):: halfgate
real,    dimension(n),intent(in ):: tdata,hdata
integer, dimension(n),intent(in ):: phof
integer,              intent(out):: doru
integer, dimension(n),intent(out):: idx,hgts
real(dp),dimension(n),intent(out):: hs
logical,              intent(out):: descending
logical,              intent(out):: FF
integer,parameter:: hour=3600 
integer          :: i,j,ii,upsign,hgs
real(dp)         :: s,gate
FF=F
if(size(hdata)/=n)stop 'In convertd; inconsistent dimensions of hdata'
if(size(tdata)/=n)stop 'In convertd; inconsistent dimensions of tdata'
if(size(hs)/=n)stop 'In convertd; inconsistent dimensions of hs'
if(size(hgts)/=n)stop 'In convertd; inconsistent dimensions of hgts'
hs=hdata
upsign=0
gate=halfgate*2
do i=1,n
   hgts(i)=2*nint(tdata(i)*hour/gate)
   if(phof(i)==5)upsign=1  
   if(phof(i)==6)upsign=-1 
enddo
doru=0
if (upsign>0) then
   doru=2
else
   doru=1
endif
if(n==1)return
if(hgts(1)>hgts(n))then 
   do i=1,n/2
      j=n+1-i
      hgs=hgts(i); hgts(i)=hgts(j); hgts(j)=hgs 
      s  =hs(i)  ; hs(i)  =hs(j)  ; hs(j)  =s   
   enddo
endif
if(upsign==1)then
   descending=F
elseif(upsign==-1)then
   descending=T
else
   descending=(hs(n)<hs(1))
   if(descending)then; upsign=-1; write(41,'("mainly DESCENDING")')
                 else; upsign=1;  write(41,'("mainly ASCENDING")')
   endif
endif
do i=1,n
   idx(i)=i
end do
do i=2,n
   do ii=1,i-1
      if (hgts(i)<hgts(ii).or. &
           (hgts(i)==hgts(ii).and.upsign*(hs(i)-hs(ii))<u0)) then  
         hgs=hgts(i);  hgts(i)=hgts(ii);   hgts(ii)=hgs   
hgts
         s  =hs(i);    hs(i)  =hs(ii);     hs(ii)  =s     
real hs
         j=idx(i);    idx(i)  =idx(ii);   idx(ii)  =j     
index idx
      end if
   end do
end do
do i=2,n
   if(hgts(i)<hgts(i-1)) then
      write(41,*)&
           'In convertd; time sequence not monotonic', i,
hgts(i),hgts(i-1)
      FF=T
      return
   end if
enddo
do i=2,n
   if(upsign*(hs(i)-hs(i-1))<u0)&
        write(41,*) 'In convertd; height sequence not monotonic'
enddo
end subroutine convertd


subroutine convertd_back(n,halfgate,wdata,tdata, &
     ws,hgts,idx,descending)
integer,              intent(in ):: n
real(dp),             intent(in ):: halfgate
integer, dimension(n),intent(in ):: hgts,idx
real(dp),dimension(n),intent(in ):: ws
logical,              intent(in ):: descending
real,dimension(n),    intent(out):: wdata
real,dimension(n),    intent(out):: tdata
integer             :: i,j,ii
real(dp)            :: s
real,   dimension(n):: wn
integer,dimension(n):: hgtn
if(size(wdata)/=n)stop 'In convertd; inconsistent dimensions of wdata'
if(size(tdata)/=n)stop 'In convertd; inconsistent dimensions of tdata'
if(size(ws)/=n)stop 'In convertd; inconsistent dimensions of ws'
if(size(hgts)/=n)stop 'In convertd; inconsistent dimensions of hgts'
do i = 1, n
   ii = idx(i); hgtn(ii) = hgts(i);  wn(ii) = ws(i)
end do
wdata=wn; tdata=hgtn*halfgate
if (descending.or.n==1) return
do i=1,n/2
   j=n+1-i
   s=tdata(i); tdata(i)=tdata(j); tdata(j)=s 
   s=wdata(i); wdata(i)=wdata(j); wdata(j)=s 
enddo
end subroutine convertd_back
