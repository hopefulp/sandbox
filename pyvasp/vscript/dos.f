!234567890
! Yong-Hyun Kim, 01May03, NREL
!     utility for projected wave functions on s- p- d- and sites
!     NEED: PROCAR, proj.in
! Yong-Hyun Kim, 09May03, NREL
!     utility for projected DOS on s- p- d- and sites
!     NEED: DOSCAR, dos.in
!
      program projected_DOS

      real, allocatable :: en(:),dos(:),tdos(:,:)
      real, allocatable :: d1(:,:),d2(:,:),d3(:,:),d4(:,:)
      integer, allocatable :: ireg(:),ion(:,:)
      real, allocatable :: rgb(:,:)

      open(unit=5,file='dos.in',form='formatted')
      read(5,*)LORBIT,ISPIN
      !if(ISPIN.ne.1.and.ISPIN.ne.-1) ISPIN=0
      if(LORBIT.eq.0) then
        nflag=4
      else
        nflag=10
      endif
      read(5,*)iprog,iorbit
      read(5,*)nreg
      allocate(ireg(nreg))
      read(5,*)(ireg(i),i=1,nreg)
      allocate(ion(nreg,maxval(abs(ireg))),rgb(nreg,3))
      ion=0
      do i=1,nreg
        if(ireg(i).gt.0) then
          read(5,*)(ion(i,j),j=1,ireg(i))
          read(5,*)rgb(i,:)
        else
          read(5,*)n1,n2
          do j=1,abs(ireg(i))
            ion(i,j)=n1+j-1
          enddo
          read(5,*)rgb(i,:)
        endif
        write(6,10)'Region ',i,':',(ion(i,j),j=1,abs(ireg(i)))
      enddo
      ireg=abs(ireg)
 10   format(a,i2,a,5000i4)
      close(5)


      open(unit=4,file='DOSCAR',form='formatted')
      read(4,*)nion
      print*,nion
      do i=1,4
      read(4,*)
      enddo
      read(4,*)fmax,fmin,ngrid,efermi
      write(6,'(a,3f10.3)')'Fermi energy: ',efermi,fmax,fmin

      allocate(en(ngrid),dos(ngrid))
      do ig=1,ngrid
        if(ISPIN.eq.-1) then
          read(4,*)en(ig),xxx,dos(ig)
        else if(ISPIN.eq.1) then
          read(4,*)en(ig),dos(ig),xxx
        else
          read(4,*)en(ig),dos(ig)
        endif
      enddo

      allocate(tdos(ngrid,nreg))
      allocate(d1(ngrid,nion),d2(ngrid,nion),d3(ngrid,nion))
      allocate(d4(ngrid,nion))

      do i=1,nion
        read(4,*,end=100)
        do ig=1,ngrid
          if(LORBIT.eq.0) then
            !read(4,*)en(ig),d1(ig,i),d2(ig,i),d3(ig,i)
            if(ISPIN.eq.-1) then
              read(4,*)en(ig),xxx,d1(ig,i),xxx,d2(ig,i),xxx,d3(ig,i)
            else if(ISPIN.eq.1) then
              read(4,*)en(ig),d1(ig,i),xxx,d2(ig,i),xxx,d3(ig,i),xxx
            else
              read(4,*)en(ig),d1(ig,i),d2(ig,i),d3(ig,i)
            endif
          else
            if(ISPIN.eq.-1) then
              read(4,*)en(ig),xxx,d1(ig,i),xxx,px,xxx,py,xxx,pz,
     1                        xxx,dd1,xxx,dd2,xxx,dd3,xxx,dd4,xxx,dd5
            else if(ISPIN.eq.1) then
              read(4,*)en(ig),d1(ig,i),xxx,px,xxx,py,xxx,pz,xxx,
     1                        dd1,xxx,dd2,xxx,dd3,xxx,dd4,xxx,dd5,xxx
            else
              read(4,*)en(ig),d1(ig,i),px,py,pz,dd1,dd2,dd3,dd4,dd5
            endif
            if(iorbit.eq.1) then
              d2(ig,i)=px
              d3(ig,i)=dd1
            else if(iorbit.eq.2) then
              d2(ig,i)=py
              d3(ig,i)=dd2
            else if(iorbit.eq.3) then
              d2(ig,i)=pz
              d3(ig,i)=dd3
            else if(iorbit.eq.4) then
              d2(ig,i)=px+py+pz
              d3(ig,i)=dd4
            else if(iorbit.eq.5) then
              d2(ig,i)=px+py+pz
              d3(ig,i)=dd5
            else
              d2(ig,i)=px+py+pz
              d3(ig,i)=dd1+dd2+dd3+dd4+dd5
            endif
          endif
        enddo
      enddo
      close(4)

      select case(iprog)
      case(1)
        d4=d1
      case(2)
        d4=d2
      case(3)
        d4=d3
      case(4)
        d4=d1+d2+d3
      case default
        STOP 'ERROR: iprog should be 1-4'
      end select


      do i=1,nreg
      do j=1,ireg(i)
        tdos(:,i)=tdos(:,i)+d4(:,ion(i,j))
      enddo
      enddo

      open(unit=8,file='dos.dat',form='formatted')
      do i=1,ngrid
       write(8,'(6f11.4)')en(i),dos(i),tdos(i,:)
      enddo
      close(8)

      deallocate(tdos,d1,d2,d3,d4)
 100  deallocate(dos,en)

      end

