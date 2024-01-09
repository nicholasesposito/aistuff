
      if ((nlv2wrt_tot.gt.1) .and. tsplines) then
         nh = 0
         do j = 1,nlv2wrt_tot
            jj = iord(j)
            if (ibfms(drinfo_accum(3,jj)).eq.0) then
               nh = nh + 1
            end if
         end do
         nh2 = nh * 2

         halfgate=30.0

         write(41,*) 'halfgate=', halfgate

         allocate(idx(nh),pof(nh))
         allocate(tdata(nh),hdata(nh),wdata(nh))
         allocate(te(nh),hgts(nh),hs(nh),dhdt(nh))
         maxita = 0
         maxitb = 0
         maxrts = 0
         maxit  = 0

         nh = 0
         do j = 1,nlv2wrt_tot
            jj = iord(j)
            if (ibfms(drinfo_accum(3,jj)).eq.0) then
               nh = nh + 1
               tdata(nh) = drinfo_accum(3,jj) 
               hdata(nh) = zevn_accum(1,jj,1) 
               pof(nh)   = nint(acft_seq_accum(2,jj))
               write(41,*) 'tdata,hdata,pof=',nh,tdata(nh),hdata(nh),
     +          pof(nh)
            end if
         end do

         call convertd(nh,halfgate,tdata,hdata,pof,
     +        doru,idx,hgts,hs,descending,FF)


         if (FF) then 
            print*,"WARNING: tspline err in utility pspl, coming out ",
     +             "of subr. convertd - use finite difference method" 
            write(41,*)"WARNING: tspline err in utility pspl, coming ",
     +                 "out of subr. convertd - use finite difference ",
     +                 "method" 
            err_tspline = 1
            go to 666
         end if
         if (descending)then
            write(41,'('' set descending'')')
         else
            write(41,'('' set ascending'')')
         endif

         call count_gates(nh,hgts(1:nh),mh)
         m = mh*2
         allocate(hgtp(m),hp(m),qbest(m),habest(m),modebest(mh))
         call best_slalom(nh,mh,doru,hgts,hs,halfgate,bigT,hgtp,hp,
     +     qbest,habest,enbest,modebest,maxita,maxitb,maxit,maxrts,FF)
          write(41,*) 'maxita,maxitb,maxit,maxrts=',maxita,maxitb,maxit,
     +     maxrts


         if (FF) then 
            print*,"WARNING: tspline err in utility pspl, coming out ",
     +             "of subr. best_slalom - use finite difference method" 
            write(41,*)"WARNING: tspline err in utility pspl, coming ",
     +                 "out of subr. best_slalom - use finite ",
     +                 "difference method" 
            err_tspline = 1
            go to 666
         end if

         call bnewton(nh,m,bigT,halfgate,hgts,hs,hgtp,habest,
     +        qbest,te(1:nh),dhdt(1:nh),FF)


         if (FF) then 
            print*,"WARNING: tspline err in utility pspl, coming out ",
     +             "of subr. bnewton - use finite difference method" 
            write(41,*)"WARNING: tspline err in utility pspl, coming ",
     +                 "out of subr. bnewton - use finite difference ",
     +                 "method" 
            err_tspline = 1
            go to 666
         end if

         call convertd_back(nh,halfgate,wdata,tdata,dhdt,hgts,idx,
     +                      descending)
         do j = 1, nh
            write(41,*) 'hgts,hs,dhdt,wdata=', j,hgts(j),hs(j),dhdt(j),
     +       wdata(j)
         end do

         nh = 0
         do j = 1,nlv2wrt_tot
            jj = iord(j)
            if (ibfms(drinfo_accum(3,jj)).eq.0) then
               nh = nh + 1
               rate_accum(jj) = wdata(nh)
               write(41,*) 'j,z,rate=',j,zevn_accum(1,jj,1),
     +          rate_accum(jj)
            end if
         end do

 666     continue

         if(allocated(idx)) deallocate(idx)
         if(allocated(pof)) deallocate(pof)
         if(allocated(tdata)) deallocate(tdata)
         if(allocated(hdata)) deallocate(hdata)
         if(allocated(wdata)) deallocate(wdata)
         if(allocated(te)) deallocate(te)
         if(allocated(hgts)) deallocate(hgts)
         if(allocated(hs)) deallocate(hs)
         if(allocated(dhdt)) deallocate(dhdt)
         if(allocated(hgtp)) deallocate(hgtp)
         if(allocated(hp)) deallocate(hp)
         if(allocated(qbest)) deallocate(qbest)
         if(allocated(habest)) deallocate(habest)
         if(allocated(modebest)) deallocate(modebest)
      end if 

      if (((nlv2wrt_tot.gt.1) .and. (.not.tsplines)) 
     +                           .or. err_tspline>0) then
        do j = 1,nlv2wrt_tot
          jj = iord(j)
          write(41,*) 'j,ord,z,t,pof=', j, jj,zevn_accum(1,jj,1),
     +    drinfo_accum(3,jj),acft_seq_accum(1,jj),acft_seq_accum(2,jj)
        end do

        do j = 1,nlv2wrt_tot
          jj = iord(j)

          jkp = 0
          jkm = 0
          jjp1 = 0
          jjm1 = 0
          if (j .eq. nlv2wrt_tot) then
             if (ibfms(drinfo_accum(3,jj)).eq.0) then
                jjp1 = jj
                jkp = j
             end if
          else
             do jk = j+1,nlv2wrt_tot
                jjp = iord(jk)
                if (jjp > nlvinprof) cycle
                if (ibfms(drinfo_accum(3,jjp)).eq.0) then
                   jjp1 = jjp
                   jkp = jk
                   exit
                end if
             end do
          end if

          if (j .eq. 1 ) then
             if (ibfms(drinfo_accum(3,jj)).eq.0) then
                jjm1 = jj
                jkm = j
             end if
          else
             do jk = j-1,1,-1
                jjm = iord(jk)
                if (jjm > nlvinprof) cycle  
                if (ibfms(drinfo_accum(3,jjm)).eq.0) then
                   jjm1 = jjm
                   jkm = jk
                   exit
                end if
             end do
          end if

          if ((jjp1 .ne. 0) .and. (jjm1 .ne. 0)) then
             dt = (drinfo_accum(3,jjp1) - drinfo_accum(3,jjm1))*3600. 

             c1_jk = 0
             c2_jk = 0
             do while ((abs(dt)<60.) .and. ((jkp+c1_jk<=nlv2wrt_tot)
     +                 .or. (jkm-c2_jk>=1)))
                jjp2 = 0
                jjm2 = 0
                c1_jk = c1_jk+1
                c2_jk = c2_jk+1
                dt_new = dt

                do while (jkp+c1_jk<=nlv2wrt_tot
     +                .and. iord(jkp+c1_jk)>nlvinprof)
                   c1_jk = c1_jk+1   
                end do
                if (jkp+c1_jk<=nlv2wrt_tot
     +                .and. iord(jkp+c1_jk)<=nlvinprof) then
                   jjp = iord(jkp+c1_jk)
                   if (ibfms(drinfo_accum(3,jjp)).eq.0) then
                      jjp2 = jjp
                      dt_new = (drinfo_accum(3,jjp2)
     +                       - drinfo_accum(3,jjm1))*3600.
                   end if
                end if
                if (abs(dt_new) >= 60.) then
                   if (jjp2 .ne. 0) jjp1 = jjp2
                   exit
                end if

                do while (jkm-c2_jk>=1 .and. iord(jkm-c2_jk)>nlvinprof)
                   c2_jk = c2_jk+1   
                end do
                if (jkm-c2_jk>=1 .and. iord(jkm-c2_jk)<=nlvinprof) then
                   jjm = iord(jkm-c2_jk)
                   if (ibfms(drinfo_accum(3,jjm)).eq.0) then
                      jjm2 = jjm
                      dt_new = (drinfo_accum(3,jjp1)
     +                       - drinfo_accum(3,jjm2))*3600.
                   end if
                end if
                if (abs(dt_new) >= 60.) then
                   if (jjm2 .ne. 0) jjm1 = jjm2
                   exit
                end if

                if ((jjp2 .ne. 0) .and. (jjm2 .ne. 0)) then
                   dt_new = (drinfo_accum(3,jjp2)
     +                    - drinfo_accum(3,jjm2))*3600.
                   if (abs(dt_new) >= 60.) then
                      if (jjp2 .ne. 0) jjp1 = jjp2
                      if (jjm2 .ne. 0) jjm1 = jjm2
                      exit
                   end if
                end if
             end do
             dt = (drinfo_accum(3,jjp1) - drinfo_accum(3,jjm1))*3600.

             zul = zevn_accum(1,jjp1,1)   
             zll = zevn_accum(1,jjm1,1) 

             if(abs(dt) .gt. 0.)  
     +          rate_accum(jj) = (zul - zll)/dt  
                                            
                                            

             write(41,*) ' fj,dt,rate_accum=',j,dt,rate_accum(jj)
             write(41,*) ''
          end if
        end do
      end if 
