      if (l_mandlvl .and. nlvinprof.gt.1) then

        do j = 1,nlv2wrt_tot
          jj = iord(j)

          nmNbtw = 0

          if(ibfms(drinfo_accum(1,jj)).ne.0 .and.
     +       ibfms(drinfo_accum(2,jj)).ne.0 .and.
     +       ibfms(drinfo_accum(3,jj)).ne.0) then

            nmNbtw = 1
            do k = j+1, nlv2wrt_tot
              kk = iord(k)
              if(ibfms(drinfo_accum(1,kk)).ne.0 .and.
     +           ibfms(drinfo_accum(2,kk)).ne.0 .and.
     +           ibfms(drinfo_accum(3,kk)).ne.0) then


                nmNbtw = nmNbtw + 1

              else
                exit
              endif
            enddo

            if(j.le.1) then
              print *
              print *, '### PROBLEM - j .le. 1 (= ',j,') in subr. ',
     +                 'sub2mem_mer, iord array underflow'
              print *, '              this should never happen
              print *
              call w3tage('PREPOBS_PREPACQC')
              call errexit(61)
            endif
            jjm1 = iord(j-1)
            jjpnmNbtw = iord(j+nmNbtw)
            pll = lvlsinprof(jjm1)
            pul = lvlsinprof(jjpnmNbtw)

            dtime_dlnp = (drinfo_accum(3,jjpnmNbtw) -
     +                    drinfo_accum(3,jjm1)) / alog(pul/pll)

            lat_pul = drinfo_accum(2,jjpnmNbtw)
            lon_pul = drinfo_accum(1,jjpnmNbtw)
            lat_pll = drinfo_accum(2,jjm1)
            lon_pll = drinfo_accum(1,jjm1)

            if(int(lon_pul*100.).eq.int(lon_pll*100.)) then
              dist_pul_pll = radius_e * abs(lat_pul-lat_pll) * deg2rad
            elseif(int(lat_pul*100.).eq.int(lat_pll*100.)) then
              dist_pul_pll = 2.0*radius_e*
     +         asin(min(1.0,abs(cos(lat_pul*deg2rad)*
     +         sin((lon_pul-lon_pll)*0.5*deg2rad))))
            else
              dist_pul_pll = 2.0*radius_e*
     +         asin(min(1.0,sqrt(
     +                           (sin((lat_pul-lat_pll)*0.5*deg2rad))**2
     +                               +  cos(lat_pul*deg2rad)*
     +                                  cos(lat_pll*deg2rad)*
     +                           (sin((lon_pul-lon_pll)*0.5*deg2rad))**2
     +                                                                 )
     +                                                                 )
     +                                                                 )
            endif

            if(int(drinfo_accum(3,jjpnmNbtw)*100000.).ne.
     +         int(drinfo_accum(3,jjm1)*100000.) .and.
     +        dist_pul_pll.ne.0.) then

              spd_pul_pll = dist_pul_pll /
     +                      abs((drinfo_accum(3,jjpnmNbtw) -
     +                           drinfo_accum(3,jjm1))*3600.)

              do k = 0,nmNbtw-1
                jjpk = iord(j+k)
                pml = lvlsinprof(jjpk)

                drinfo_accum(3,jjpk) = drinfo_accum(3,jjm1)
     +                                 dtime_dlnp*alog(pml/pll)

                dist2pml = spd_pul_pll *
     +                   abs(drinfo_accum(3,jjpk)-drinfo_accum(3,jjm1))*
     +                     3600.

                drinfo_accum(2,jjpk) = drinfo_accum(2,jjm1)
     +                                 dist2pml/dist_pul_pll*
     +           (drinfo_accum(2,jjpnmNbtw)-drinfo_accum(2,jjm1))

                drinfo_accum(1,jjpk) = drinfo_accum(1,jjm1)
     +                                 dist2pml/dist_pul_pll*
     +           (drinfo_accum(1,jjpnmNbtw)-drinfo_accum(1,jjm1))

              enddo
            else



              delx = (drinfo_accum(1,jjpnmNbtw) -
     +                drinfo_accum(1,jjm1))/(nmNbtw+1)
              dely = (drinfo_accum(2,jjpnmNbtw) -
     +                drinfo_accum(2,jjm1))/(nmNbtw+1)

              do k = 0,nmNbtw-1
                jjpk = iord(j+k)
                pml = lvlsinprof(jjpk)
                drinfo_accum(1,jjpk) =
     +           drinfo_accum(1,jjm1) + (k+1)*delx
                drinfo_accum(2,jjpk) =
     +           drinfo_accum(2,jjm1) + (k+1)*dely
                drinfo_accum(3,jjpk) = drinfo_accum(3,jjm1)
     +                                 dtime_dlnp*alog(pml/pll)

              enddo
            endif
          endif

        enddo
      endif


      jjmaxp = iord(nlv2wrt_tot)
      jjminp = iord(1)
      if(nlv2wrt_tot.eq.1) then
        hdr2wrt(6) = 300 + mod(int(hdr2wrt(6)),100)   
                                                      
                                                      
      elseif(nlv2wrt_tot.gt.1 .and.
     +      (c_qc_accum(jjmaxp)(11:11).eq.'a' .or.
     +       c_qc_accum(jjmaxp)(11:11).eq.'A')) then  
        hdr2wrt(6) = 400 + mod(int(hdr2wrt(6)),100)   
                                                      
                                                      
                                                      

        hdr2wrt(2) = drinfo_accum(1,jjmaxp)
        hdr2wrt(3) = drinfo_accum(2,jjmaxp)
        hdr2wrt(4) = drinfo_accum(3,jjmaxp)
        hdr2wrt(5) = elv_accum(1,jjmaxp)
        hdr2wrt(12) = rpt_accum(1,jjmaxp)
        hdr2wrt(13) = tcor_accum(1,jjmaxp)

      elseif(nlv2wrt_tot.gt.1 .and.
     +      (c_qc_accum(jjmaxp)(11:11).eq.'d' .or.
     +       c_qc_accum(jjmaxp)(11:11).eq.'D')) then  
        hdr2wrt(6) = 500 + mod(int(hdr2wrt(6)),100)   
                                                      
                                                      
                                                      

        hdr2wrt(2) = drinfo_accum(1,jjminp)
        hdr2wrt(3) = drinfo_accum(2,jjminp)
        hdr2wrt(4) = drinfo_accum(3,jjminp)
        hdr2wrt(5) = elv_accum(1,jjminp)
        hdr2wrt(12) = rpt_accum(1,jjminp)
        hdr2wrt(13) = tcor_accum(1,jjminp)

      endif

      hdr2wrt(10) = bmiss
      hdr2wrt(11) = bmiss

