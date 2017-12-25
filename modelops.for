c-----------------------------------------------------------------------
c
c   Apply algebriac operations on multiple datafiles and data columns 
c   for files in GSLIB Geo-EAS format.  
c   Specific any number of data files and columns as inputs they are 
c   stored internally as data columns.  Apply any number of algebriac
c   operations resulting each in a new data column.  The select an
c   output file with any data column.  The program allows for fast and 
c   flexible manipulation of models that are too large to be treated
c   with Excel, Python or R efficiently.
c
c   By: Michael Pyrcz, University of Texas at Austin, @GeostatGuy
c-----------------------------------------------------------------------
      use       msflib
      parameter(MAXLEN=512,EPSLON=0.00001,VERSION=1.000)

      real      var(500),dataloc(2,20),tmin,tmax
	integer   test
      character outfl*512,str*512,strlin*132
	logical   testfl
      data      lin/1/,lout/2/

	character, allocatable :: datafl(:)*512
      integer, allocatable   :: opmat(:,:),outlist(:),maskmat(:)
	real, allocatable      :: datamat(:,:)

c
c Note VERSION number before anything else:
c
      write(*,9999) VERSION
 9999 format(/' modelops Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'modelops.par         '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'modelops.par           ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a40)',end=98) str
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c

      read(lin,*,err=98) nfile,nop,tmin,tmax
      write(*,*) ' nfile,nop,tmin,tmax = ',nfile,nop,tmin,tmax

      read(lin,*,err=98) icalc
      write(*,*) ' icalc = ',icalc

c
      allocate(datafl(nfile),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(opmat(4,nop),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c

      idata = 1
      do ifile = 1, nfile
       read(lin,'(a512)',err=98) datafl(ifile)
       call chknam(datafl(ifile),512)
       write(*,'(a40)') ' Input file - ',datafl(ifile)
	 read(lin,*,err=98) ndata
	 backspace(lin)
	 read(lin,*,err=98) ndata, (dataloc(2,i), i=idata,idata+ndata-1)
       write(*,*) 'ndata,cols',
     +  ndata,(dataloc(2,i),i=idata,idata+ndata-1)
	 do i = idata,idata+ndata-1
	  dataloc(1,i) = ifile
	 end do
	 idata = idata + ndata 
      end do
      ndata = idata-1

      read(lin,*,err=98) nx,ny,nz
      write(*,*) ' nx,ny,nz = ',nx,ny,nz

      do iop = 1, nop
	 read(lin,*,err=98) (opmat(i,iop),i=1,4)
	 write(*,*) ' d1,d2,op,d3 = ',(opmat(i,iop),i=1,4)
      end do

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) osx,oex
      write(*,*) 'Output: ix=osx,..,oex = ',osx,oex

      read(lin,*,err=98) osy,oey 
      write(*,*) 'Output: iy=osy,..,oey = ',osy,oey

      read(lin,*,err=98) osz,oez
      write(*,*) 'Output: iz=osz,..,oez = ',osz,oez

      read(lin,*,err=98) nout
c
      allocate(outlist(nout),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c

      backspace(lin)
      read(lin,*,err=98) nout,(outlist(i),i=1,nout)
      write(*,*) ' output columns = ',nout,(outlist(i),i=1,nout)

      close(lin)

c
c Set some parameters
c

      ncell = nx*ny*nz
      nval = 1
      do iop = 1, nop
        nval = max(nval,opmat(4,iop))
      end do
	 

c
      allocate(datamat(nx,ny,nz,nval),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c
      allocate(maskmat(nx,ny,nz),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed',
     +                  ' due to insufficient memory.'
                  stop
            end if
c

c
c Open the input files and read them in
c
      datamat=0.0
      do idata = 1, ndata
        write(*,*) 'Reading in the Data #',idata
	inquire(file=datafl(dataloc(1,idata)),exist=testfl)
        if(.not.testfl) then
            write(*,*) 'Data file ',datafl(dataloc(1,idata)),
     +	  ' does not exist! - assuming constant values from icol'
          do iz = 1, nz
	    do iy = 1, ny
	      do ix = 1, nx
	       datamat(ix,iy,iz,idata) = dataloc(2,idata)
	      end do
            end do
          end do       	 
	else
          open(lin,file=datafl(int(dataloc(1,idata))),status='OLD')
          read(lin,*,err=98)
          read(lin,*,err=98) nvari
          do i=1,nvari
            read(lin,*,err=98)
          end do
          do iz = 1, nz
	    do iy = 1, ny
	      do ix = 1, nx
	        read(lin,*,err=98) (var(i),i=1,nvari)
	        val=var(int(dataloc(2,idata)))
	        if(val.lt.tmin.or.val.gt.tmax) then
	          maskmat(icell) = 1
                end if
	        datamat(ix,iy,iz,idata) = val
	      end do
            end do
	  end do
	  close(lin)
        end if
      end do

c
c Complete the modelops operations
c

      do iop = 1, nop
        write(*,*) 'Working on operation #',iop
        d1 = opmat(1,iop)
	d2 = opmat(2,iop)
	op = opmat(3,iop)
	d3 = opmat(4,iop)
	if(op.eq.1) then
	 do iz = 1, nz
	  do iy = 1, ny
	   do ix = 1, nx
	    datamat(ix,iy,iz,d3) = datamat(ix,iy,iz,d1)+
     +       datamat(ix,iy,iz.d2)
	   end do
          end do
         end do
	else if(op.eq.2) then
	 do iz = 1, nz
	  do iy = 1, ny
	   do ix = 1, nx
	    datamat(ix,iy,iz,d3) = datamat(ix,iy,iz,d1)-
     +       datamat(ix,iy,iz,d2)
	   end do
          end do
         end do
	else if(op.eq.3) then
	 do iz = 1, nz
	  do iy = 1, ny
	   do ix = 1, nx
	    datamat(ix,iy,iz,d3) = datamat(ix,iy,iz,d1)*
     +       datamat(ix,iy,iz,d2)
	   end do
          end do
         end do 
	else if(op.eq.4) then
	 do iz = 1, nz
	  do iy = 1, ny
	   do ix = 1, nx
	    if(datamat(ix,iy,iz,d2).ne.0.0) then
	     datamat(ix,iy,iz,d3) = datamat(ix,iy,iz,d1)/
     +        datamat(ix,iy,iz,d2)
	    else
	     datamat(ix,iy,iz,d3) = -99999.9
	    end if
	   end do
          end do
         end do 	
	else if(op.eq.5) then
	 do iz = 1, nz
	  do iy = 1, ny
	   do ix = 1, nx
	    datamat(d3,icell) = datamat(d1,icell)** 
     +       datamat(d2,icell)
	   end do
          end do
         end do
        end if
c        else if(op.eq.6) then
      end do

c
c Write out the output
c
      open(lout,file=outfl,status='UNKNOWN')
      write(lout,*) 'modelops_Output'
      write(lout,*) nout
      do iout = 1, nout
	write(lout,*) 'Data',iout
      end do

      do icell = 1, ncell
        if(maskmat(icell).eq.1) then
	  do iout = 1, nout
	    datamat(outlist(iout),icell) = tmin-100.0
	  end do
	end if
	write(lout,531) (datamat(outlist(iout),icell),iout=1,nout)
      end do
      close(lout)

531   format(30(f15.5,1x))

c
c Finished:
c
      write(*,9998) VERSION
 9998 format(/' modelops Version: ',f5.3, ' Finished'/)
      stop
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end

      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='modelops.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for modelops',/,
     +       '                  ***********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('2 2 -1.0 100.0                ',
     +       '-nfile, nop, tmin, tmax')
      write(lun,12)
 12   format('data1.dat                     ',
     +       '-datafl#1')
      write(lun,13)
 13   format('1 1                           ',
     +       '-datafl#1: ndata, cols')
      write(lun,14)
 14   format('data2.dat                     ',
     +       '-datafl #2')
      write(lun,15)
 15   format('1 1                           ',
     +       '-datafl#1: ndata, cols')
       write(lun,16)
 16   format('10 10 1                       ',
     +       '-nx,ny,nz')
      write(lun,17)
 17   format('1 2 1 3                       ',
     +       '-operation#1: d1, d2, op, d3')
      write(lun,18)
 18   format('1 3 1 4                       ',
     +       '-operation#2: d1, d2, op, d3')
      write(lun,19)
 19   format('modelops.out                  ',
     +       '-output file')
      write(lun,20)
 20   format('4 1 2 3 4                     ',
     +       '-nout, values')
      write(lun,*)
      write(lun,21)
 21   format('1:d1+d2,2:d1-d2,3:d3*d4,4:d3/d4,5:d3**d4          ')

      close(lun)
      return
      end

