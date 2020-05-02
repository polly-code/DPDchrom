implicit none
include "mpif.h"

common /cell/ dlx,dly,dlz,rcut,kpbc !box size on each node (dlx,dly,dlz), cut radius (rcut), type of PBC(kpbc)
common /cella/ dlxa,dlya,dlza !complete cox size
common /sizes/ sz !size of main arrays (like r and v)
common /bxsize/ ilx,ily,ilz !number of box bins 
common /na/ natms,natmsfan !current number of "real" particles(natms), current number of right "ghost" particles
common /nab/ natmsfanbond !current number of left "ghost" particles. necessary for bond forces calculation
common /datatech/ alpha,sigma,gamma !dpd parameters 
common /time/ dt ! timestep
common /comm/ ssz,fsz !size of communication arrays, ssz is for moving atoms between nodes, fsz - for "ghost" atoms 
common /npall/ npall !total number of atoms in the system
common /shwi/ shwi !transfer radius for "ghost" atoms 
common /correl/ corr,bcorr,bcorrt !arrays that store global atom numbers. corr converts atom number on the node to the global number, bcorr does the opposite, bcorrt shows atom state ("real,"ghost" or on the other node)
common /den/ rho !density 
common /bondf/ k_eq,k_bond !bond parameters - eq.length, stiffness
common /fname/ msi !name of input file
common /mpidpd/ rank, nproc ! node number and total number of nodes
common /rxyz/ rx,ry,rz !coordinates
common /fxyz/ fx,fy,fz !forces
common /vxyz/ vx,vy,vz !velocities
common /vwxyz/ vwx,vwy,vwz !waved velocities(for integration) 
common /ptype/ kindp !current particle types
common /sndrcv/ sndr,sndl,fanl,fanr,rcv,rcvf !communication arrays. sndr,sndl - for moving particles, fanl,fanr - for "ghost", rcv,rcvf - for receive 
common /funnumb/ npl,fann !npl - number of "ghost" atoms than have been sent to left node, fann - numbers of these particles. necessary for receiving forces actiong on atoms close to the node boundary from left node
common /commf/ fsend,rcou !array that stores forces to be sent to the right node (see above); rcou - detailed number of ghost particles
common /bnd/ cn !main bond array
common /rst/ nbond !number of bonds
common /typel/ typelist,ntype !type of the atoms that present in the system and number of such types
common /ptypeall/ kindpall !particle types
common /top/ cln,nrow,ncol ! number of nodes in 1 row,row and column number of node
common /prs/ press ! pressure
common /vale/ val !valencies
common /bmult/ mult !bond snd rcv arrays size multiplier
common /lambd/ lambda
common /bio_input/ resolution, geo_file, contacts ! input parameters
common /reinit/ rein ! reinitialize variable for speeding up the calculations

real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:) !size ssz,4;ssz,4;fsz,4;fsz,4;ssz*8,fsz*8
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) ! size sz,3
real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) ! size sz,2
real*4,pointer :: fx(:),fy(:),fz(:) ! size sz
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
integer*1,pointer :: kindp(:,:) ! size sz,3
integer*4,pointer :: corr(:,:),bcorr(:)! size sz,3 and npall
integer*1,pointer :: bcorrt(:) !size npall
integer*4,pointer :: fann(:,:) !size fsz/8,4
integer*4,pointer :: contacts(:,:) !size nbnd,2
real*4,pointer :: fsend(:) ! size 3*fsz/2
integer*4,pointer :: cn(:,:) !size npall,0:maximum number of bonds
integer*1,pointer :: kindpall(:) ! size npall
integer*1,pointer :: val(:) !size npall
real*4,pointer :: alpha(:,:) !size 10,10
real*4,pointer :: k_eq(:,:),k_bond(:,:) !size 10,10

integer*4 nbond
integer*4 rank,nproc,cln,nrow,ncol
real*4 uni
integer*4 sz,ssz,fsz,npall
integer*4 ilx,ily,ilz
integer ierr
integer*4 steps1,steps2,qstp,nstp
integer*4 kpbc,rst
real*4 dlx,dly,dlz,rcut,dlxtar
real*4 dlxa,dlya,dlza
integer*4 rcou(4)
integer*4 natms,natmsfan,npl(4),natmsfanbond,mult
real*4 sigma,dt,gamma
real*4 shwi
real*4 rho
character*8  msi 
real*4 vall(3),vl(1)
real*8 start
integer*1 typelist(10),ntype
integer*4 i,j,curdef
real*4 press(3)
real*4 lambda
integer*8 resolution
character*256 geo_file
integer*4 rein

! MPI initialization
call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)

! randon number generator initialization
vl=uni()
call read_from_command_line()

call barrier

! read input script
!call rconf(steps1,steps2,steps3,qstp,rst,rststp,stepsdef,ndef,tau,tai,def_res)
call rconf(steps1,steps2,qstp)
! read restart.dat file
call readrst()

! allocate arrays
allocate (vwx(sz,2),vwy(sz,2),vwz(sz,2),fx(sz),fy(sz),fz(sz))
allocate (fann(fsz/8,4))
allocate (fsend(3*fsz/2))
allocate (fanl(max(fsz,ssz),4),fanr(max(fsz,ssz),4),rcvf(max(fsz,ssz)*8))
sndr => fanr
sndl => fanl
rcv => rcvf

! generate initial ghost beads arrays
call initial()

! output information
if(rank==0) then
    write(*,*)
    write(*,*) '*********************************************'
    write(*,*)
    write(*,*) '             DPDchrom'
    write(*,*)
    write(*,*) 'Some technical information'
    write(*,*) 'Number of nodes              :',nproc
    write(*,*) 'Number of nodes in 1 row     :',cln
    write(*,*)
    write(*,*) 'System information'
    write(*,*) 'Box size X Y Z               :',dlxa,dlya,dlza
    write(*,*) 'Box size at each node        :',dlx,dly,dlz
    write(*,*) 'DPD density                  :',rho  
    write(*,*) 'Total number of particles    :',npall
    write(*,*) 'Types of particles           :',typelist(1:ntype)
    write(*,*)
    write(*,*) 'Integration information'
    write(*,*) 'Timestep                     :',dt
    write(*,*) 'Number of steps in 1st stage :',steps1
    write(*,*) 'Number of steps in 2nd stage :',steps2
    write(*,*)
    write(*,*) 'DPD Parameters'
    write(*,*) 'Sigma                        :',sigma
    write(*,*) 'Gamma                        :',gamma
    do i=1,ntype
        do j=i,ntype
            write(*,*) 'DPD coeff btw',typelist(i),'and',typelist(j),'is',alpha(typelist(i),typelist(j))
        end do
    end do

    write(*,*)
    write(*,*) '*********************************************'
    write(*,*)
end if

! calculate forces
fsend=0
call forces()
call forcesbond()

if (rank==0) start=mpi_wtime()

! start 1st stage with soft potentials
alpha=alpha/5.0
k_bond=k_bond/5.0
dt=dt/10 !timestep for equilibration
rein=1 ! reinit all variables with dt
do nstp=1,steps1 ! stage 1 - equilibration 
    
    ! integrate beads positions and velocities
    call intr()

    ! calculate forces
    call forces()
    call forcesbond()

    ! integrate beads velocoties
    call intv()
    
    ! control output
    if (mod(nstp,qstp)==0) then
        vl=(sum(vz*vz)+sum(vy*vy)+sum(vx*vx))
        vl=vl/natms/3
        call gsum(vl,1)
        vl=vl/nproc
        vall(1)=sum(vx)
        vall(2)=sum(vy)
        vall(3)=sum(vz)
        call gsum(vall,3)
        call pressure()
        if (rank==0) write(*,*)'stage 1', nstp,vl,maxval(abs(vall)),sum(press)/3.0
    end if
    
    ! control  CM velocily
    if (mod(nstp,100)==0) call vcontrol()
end do

! finish 1st pahse
alpha=alpha*5.0
k_bond=k_bond*5.0
dt=dt*10
rein=1 ! reinit all variables with dt
! output time
if (rank==0) write(*,*) 'stage 1 took',mpi_wtime()-start,'s'
! start 2nd stage - calculation
if (rank==0) start=mpi_wtime()
do nstp=1,steps2

    ! integrate beads positions and velocities
    call intr()

    ! calculate forces
    call forces()
    call forcesbond()

    ! integrate beads velocoties
    call intv()

    ! control output
    if (mod(nstp,qstp)==0) then
        vl=(sum(vz*vz)+sum(vy*vy)+sum(vx*vx))
        vl=vl/natms/3
        call gsum(vl,1)
        vl=vl/nproc
        vall(1)=sum(vx)
        vall(2)=sum(vy)
        vall(3)=sum(vz)
        call gsum(vall,3)
        call pressure()
        if (rank==0) write(*,*)'stage 2', nstp,vl,maxval(abs(vall)),sum(press)/3.0
    end if

    if (mod(nstp,10000)==0) call writerst()
    ! control  CM velocily
    if (mod(nstp,100)==0) call vcontrol()
end do

! finish 2nd stage 
! output time
if (rank==0) write(*,*) 'stage 2 took',mpi_wtime()-start,'s'

! output restart file
call writerst()

! delete arrays
call barrier
nullify(rx,ry,rz,vwx,vwy,vwz,fx,fy,fz,vx,vy,vz,kindp,corr,bcorr,bcorrt,sndr,sndl,fanl,rcv,rcvf,fann,fsend)

! end mpi
call mpi_finalize(ierr)

end

subroutine read_from_command_line()
! #####################################################################
! #                                                                   #
! #   subroutine -1:                                                  #
! #   read geo file                                                   #
! #                                                                   #
! #####################################################################
implicit none
save
common /bio_input/ resolution, geo_file, contacts ! input parameters

character*256 geo_file, resstr
integer*8 resolution, num_of_args
integer*4,pointer :: contacts(:,:)

num_of_args = IARGC()
if (num_of_args .ne. 2) stop 'Number of arguments not equals 2! The first one is filename, the second one is resolution. For example ". ./dpd name.geo 100000"'
CALL GETARG (1,geo_file)
CALL GETARG (2,resstr)
read(resstr,'(I18)') resolution

end

subroutine read_geo()
! #####################################################################
! #                                                                   #
! #   subroutine 0:                                                   #
! #   read geo file                                                   #
! #                                                                   #
! #####################################################################
implicit none
save
common /bio_input/ resolution, geo_file, contacts ! input parameters
common /cella/ dlxa,dlya,dlza
common /den/ rho

integer*8 :: maxChrLength(100), minChrLength(100)
integer*4,pointer :: contacts(:,:), preset(:)
integer*4,pointer :: cn2(:,:)
integer*4, pointer :: arrLengthsOfChains(:)
real*4,pointer :: rxt(:),ryt(:),rzt(:)
character*32 :: chrName(100)

integer*8 resolution
integer*8 pos1,pos2
integer*4 nlines
integer*4 numberOfChr, iter
integer*4 zero, one, two
integer*4 i,j, i1, i2,nbmax2,natmsall,nbondsall
integer*4 totalNumberOfSolvent,totalNumberOfChainBeads
integer*4 globalPos1, globalPos2, reason, iter2, empt_var
real*4 dlxa,dlya,dlza, rho
real*4 rndposx, rndposy, rndposz
real*4 uni
character*32 chr1,chr2
character*32 strand1,strand2
character*256 geo_file, trash_line, test_line(6)
logical firstFlag, secondFlag, skip_1st_line

data firstFlag/.false./
data secondFlag/.false./
data skip_1st_line/.true./

zero=0
one=1
two=2
rho=3.0
numberOfChr=0
nlines=-1 !first line is a header
chrName='empty'

open(14, file = geo_file, status = 'old')

do
    if (skip_1st_line) then
        skip_1st_line=.false.
        read(14,'(a)') trash_line
        cycle
    end if
    read(14,*,iostat=reason)chr1,chr2,pos1,pos2,strand1,strand2
    if (reason>0) then
        stop 'something wrong with the file geo'
    else if (reason<0) then !end of file
        exit
    end if

    do i=1,100
        if ((.not. firstFlag) .and. (chrName(i)=='empty')) then
            chrName(i)=chr1
            maxChrLength(i)=pos1
            minChrLength(i)=pos1
            numberOfChr=numberOfChr+1
            firstFlag=.true.
        else if (.not.firstFlag .and.chrName(i)==chr1) then
            if (pos1>maxChrLength(i))maxChrLength(i)=pos1
            if (pos1<minChrLength(i))minChrLength(i)=pos1
            firstFlag=.true.
        end if
        if (.not.secondFlag .and. chrName(i)=='empty') then
            chrName(i)=chr2
            maxChrLength(i)=pos2
            minChrLength(i)=pos2
            numberOfChr=numberOfChr+1
            secondFlag=.true.
        else if (.not.secondFlag .and.chrName(i)==chr2) then
            if (pos2>maxChrLength(i))maxChrLength(i)=pos2
            if (pos2<minChrLength(i))minChrLength(i)=pos2
            secondFlag=.true.
        end if
        if (firstFlag .and. secondFlag) then
            firstFlag=.false.
            secondFlag=.false.
            exit
        end if
    end do
    nlines=nlines+1
end do
close(14, status = 'keep')
open(14, file = geo_file, status = 'old')
allocate(arrLengthsOfChains(numberOfChr))
do i=1,numberOfChr
    arrLengthsOfChains(i)=CEILING(real(maxChrLength(i)-minChrLength(i))/real(resolution))
end do
totalNumberOfChainBeads = sum(arrLengthsOfChains)
allocate(contacts(totalNumberOfChainBeads**2,2))
contacts=0
iter=1
iter2=1

do i=1,numberOfChr
write(*,*)'Length of chromosome',i,arrLengthsOfChains(i)
    do j=1,arrLengthsOfChains(i)-1
        contacts(iter,1)=iter2
        contacts(iter,2)=iter2+1
        iter=iter+1
        iter2=iter2+1
        if (j==arrLengthsOfChains(i)-1) iter2=iter2+1
    end do
end do

skip_1st_line=.true.
iter2=iter-1
allocate(preset(numberOfChr))

preset(1)=0
do i=2,numberOfChr
    preset(i)=sum(arrLengthsOfChains(1: i-1))
end do

do
    if (skip_1st_line) then
        skip_1st_line=.false.
        read(14,'(a)') trash_line
        cycle
    end if
    read(14,*,iostat=reason)chr1,chr2,pos1,pos2,strand1,strand2
    if (reason>0) then
        stop 'something wrong with the file geo'
    else if (reason<0) then !end of file
        exit
    end if
    firstFlag=.true.
    do i=1,numberOfChr
        if (chr1==chrName(i)) globalPos1=preset(i)+ceiling(real(pos1-minChrLength(i)+1)/real(resolution))
        if (chr2==chrName(i)) globalPos2=preset(i)+ceiling(real(pos2-minChrLength(i)+1)/real(resolution))
    end do
    if (abs(globalPos1-globalPos2)>1) then
        do j=iter2,iter
            if ((contacts(j,1)==globalPos1 .and. contacts(j,2)==globalPos2).or.(contacts(j,1)==globalPos2 .and. contacts(j,2)==globalPos1)) then
                firstFlag=.false.
                exit
            end if
        end do
        if (firstFlag) then
            contacts(iter,1)=globalPos1
            contacts(iter,2)=globalPos2
            iter=iter+1
        end if
        firstFlag=.true.
    end if
end do
close(14, status = 'keep')
!calculate the size of simulation box
dlxa=((totalNumberOfChainBeads/rho)**(1./3.))*2+1
dlya=dlxa
dlza=dlxa
natmsall=dlxa**3*rho
totalNumberOfSolvent=natmsall-totalNumberOfChainBeads

allocate(cn2(natmsall,0:0)) 
cn2=0
nbondsall=iter-1
do i = 1, nbondsall

    cn2(contacts(i,1),0)=cn2(contacts(i,1),0)+1
    cn2(contacts(i,2),0)=cn2(contacts(i,2),0)+1
end do

nbmax2=maxval(cn2(:,0))
nullify (cn2)
allocate (cn2(natmsall,0:nbmax2))

! constuct connection matrices for all bonds
cn2=0
do i = 1, nbondsall
    i1=contacts(i,1)
    i2=contacts(i,2)
    cn2(i1,0) = cn2(i1,0) + 1
    cn2(i1,cn2(i1,0)) = i2
    cn2(i2,0) = cn2(i2,0) + 1
    cn2(i2,cn2(i2,0)) = i1
end do
write(*,*) 'start to write restart, atoms, bonds:', natmsall, nbondsall
!write the data to restart file
open(3, file = 'restart.dat',status='replace' )
write(3,'(i8,f8.4)') natmsall, rho
write(3,'(3f20.12)') dlxa,dlya,dlza

iter=1

do i=1,numberOfChr
    allocate(rxt(arrLengthsOfChains(i)),ryt(arrLengthsOfChains(i)),rzt(arrLengthsOfChains(i)))
    call rndplace( arrLengthsOfChains(i),rxt,ryt,rzt,cn2,(dlxa-1)/2,(dlya-1)/2,(dlza-1)/2, 2 )
    rndposx = 0
    rndposy = 0
    rndposz = 0
    do j=1,arrLengthsOfChains(i)
        write(3,'(i8,2i4,3f14.6)')iter, zero, one, rxt(j)+3, ryt(j)+3, rzt(j)+3
        iter=iter+1
    end do
    nullify(rxt,ryt,rzt)
end do

do i=1,totalNumberOfSolvent
    write(3,'(i8,2i4,3f14.6)')iter, zero, two, dlxa*(uni()), dlya*(uni()), dlza*(uni())
    iter=iter+1
end do

! write information about bonds
write(3,'(a9,i9)') ' bonds:  ',nbondsall
do i=1,nbondsall
    write(3,'(2i8)') contacts(i,1), contacts(i,2)
end do
write(3,'(a9,i9)') ' angles:  ',zero
close(3)
end

 subroutine forces()
! #####################################################################
! #                                                                   #
! #   subroutine 1:                                                   #
! #   calculate all forces: conservative, dissipative, and random     #
! #                                                                   #
! #####################################################################
implicit none
save
      
common /cell/ dlx,dly,dlz,rcut,kpbc
common /sizes/ sz
common /bxsize/ ilx,ily,ilz
common /na/ natms,natmsfan
common /rxyz/ rx,ry,rz
common /fxyz/ fx,fy,fz
common /vwxyz/ vwx,vwy,vwz
common /ptype/ kindp
common /shwi/ shwi 
common /prs/ press 

integer*4 sz
integer*4 ilx,ily,ilz        
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 
real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) 
real*4,pointer :: fx(:),fy(:),fz(:) 
integer*1,pointer :: kindp(:,:)
integer*4,pointer :: fmap(:,:,:) 
real*4,pointer :: fmapc(:,:) 

logical init
integer*4 link(sz,2),lct(ilx+1,0:ily+1,ilz)
integer*4 natms,natmsfan
integer*4 ncells,ix,iy,iz,jx,jy,jz
real*4 dlx,dly,dlz,rcut,rcsq,shwi
integer*4 i,j,ic,kc,ii
integer*4 kpbc
real*4 xd,yd,zd,rsq,cz
integer*4 arrtype
real*4 press(3)

data init/.true./

! set initial values
if (init) then
    init=.false.
    rcsq=rcut**2 
    allocate (fmap(3,18,(ilx+1)*(ily+1)*(ilz+1)),fmapc(18,(ilx+1)*(ily+1)*(ilz+1)))
end if

! create subcells map        
ncells=ilx*ily*ilz
call fmaps(fmap,fmapc)

! zero arrays
fx=0
fy=0
fz=0
press=0
link=0
lct=0

! create linked list for real particles
do i=1,natms
    ix=min(int(rx(i,1)/rcut),ilx-1)+1
    iy=min(int(ry(i,1)/rcut),ily-1)+1
    iz=min(int(rz(i,1)/rcut),ilz-1)+1

    j=lct(ix,iy,iz)
    lct(ix,iy,iz)=i
    link(i,1)=j
end do
 
! create linked list for ghost particles        
do i=1,natmsfan

    if (rx(i,2)>=(dlx+rcut)) then
        cycle
    elseif (rx(i,2)<dlx) then
        ix=min(int(rx(i,2)/rcut),ilx-1)+1
    elseif (rx(i,2)<(dlx+rcut)) then
        ix=ilx+1
    end if
          
    if (ry(i,2)>=(dly+rcut).or.ry(i,2)<=(-rcut)) then
        cycle
    elseif (ry(i,2)<0) then
        iy=0  
    elseif (ry(i,2)<dly) then
        iy=min(int(ry(i,2)/rcut),ily-1)+1
    elseif (ry(i,2)<(dly+rcut)) then
        iy=ily+1
    end if
          
    iz=min(int(rz(i,2)/rcut),ilz-1)+1

    j=lct(ix,iy,iz)
    lct(ix,iy,iz)=i
    link(i,2)=j
end do

! primary loop over subcells
iz=1
iy=1
ix=1
        
do ic=1,ncells
    ii=lct(ix,iy,iz)
    if (ii>0) then
           
        ! secondary loop over subcells
        do kc=1,17
            i=ii
            cz=fmapc(kc,ic)
            jx=fmap(1,kc,ic)
            jy=fmap(2,kc,ic)
            jz=fmap(3,kc,ic)
            if (jx<0) exit
            j=lct(jx,jy,jz)

            ! ignore if empty
            if (j>0) then
                arrtype=1
                if (jx>ilx.or.jy>ily) arrtype=2
                do while(i/=0)
                    ! check if the subcells are identical 
                    if ((jx==ix).and.(jy==iy).and.(jz==iz)) j=link(i,1)
                    if (j>0) then
                        do while(j/=0)

                            ! distance in space
                            xd=-(rx(j,arrtype)-rx(i,1))
                            yd=-(ry(j,arrtype)-ry(i,1))
                            zd=-(rz(j,arrtype)-rz(i,1)+cz)
                            rsq = xd*xd + yd*yd + zd*zd

                            ! test the distance
                            if (rcsq>rsq) call fcdrij(i,j,rsq,xd,yd,zd,arrtype)

                            ! go to the next particle (second subcell)
                            j=link(j,arrtype)
                        end do
                    end if

                    ! return to the first particle (second subcell)
                    j=lct(jx,jy,jz)

                    ! go to the next particle (first subcell)
                    i=link(i,1)
                end do
            end if
        end do
    end if

    ! move to the next subcell 
    ix=ix+1
    if (ix>ilx) then
        ix=1
        iy=iy+1
        if(iy>ily)then
            iy=1
            iz=iz+1
        end if
    end if
end do

! send forces acting on right "ghost" particles
call fsendunpck()

end
      
subroutine fmaps(fmap,fmapc)
! #####################################################################
! #                                                                   #
! #   subroutine 2:                                                   #
! #   create map of cubcells neighbours                               #
! #   for forces                                                      #
! #                                                                   #
! #####################################################################
implicit none

common /bxsize/ ilx,ily,ilz
common /cell/ dlx,dly,dlz,rcut,kpbc

real*4 dlx,dly,dlz,rcut
integer*4 kpbc
integer*4 ilx,ily,ilz
integer*4 ic,kc
integer*4 nsbcll,ncells,ix,iy,iz,jx,jy,jz
real*4 cy,cz
integer*4 fmap(3,18,*)
real*4 fmapc(18,*)

integer*4 nix,niy,niz,nix1,niy1,niz1,nix2,niy2,niz2,nix3,niy3,niz3,nix4,niy4,niz4
dimension nix(18),niy(18),niz(18),nix1(14),niy1(14),niz1(14),nix2(17),niy2(17),niz2(17),nix3(11),niy3(11),niz3(11),nix4(14),niy4(14),niz4(14)

data nix1/0,0,0,1, 0,0, 1, 1,1,1, 1, 1, 1,1/
data niy1/0,0,1,0, 1,1,-1, 0,0,1,-1,-1, 1,1/
data niz1/0,1,0,0,-1,1, 0,-1,1,0,-1, 1,-1,1/

data nix2/0,0,0,1, 0,0, 1, 1,1,1, 1, 1, 1,1,-1,-1,-1/
data niy2/0,0,1,0, 1,1,-1, 0,0,1,-1,-1, 1,1, 1, 1, 1/
data niz2/0,1,0,0,-1,1, 0,-1,1,0,-1, 1,-1,1, 1, 0,-1/

data nix3/0,0,0,1, 0,0, 1,1,1, 1,1/
data niy3/0,0,1,0, 1,1, 0,0,1, 1,1/
data niz3/0,1,0,0,-1,1,-1,1,0,-1,1/

data nix4/0,0,0,1, 0,0, 1,1,1, 1,1,-1,-1,-1/
data niy4/0,0,1,0, 1,1, 0,0,1, 1,1, 1, 1, 1/
data niz4/0,1,0,0,-1,1,-1,1,0,-1,1, 1, 0,-1/

iz=1
iy=1
ix=1
fmap(:,:,1:ilx*ily*ilz)=-1
fmapc(:,1:ilx*ily*ilz)=0
ncells=ilx*ily*ilz

! loop over subcells
do ic=1,ncells

    ! case 1 - thin domain along X
    if (ilx==1) then
        nix(1:14)=nix1
        niy(1:14)=niy1
        niz(1:14)=niz1
        nsbcll=14
    goto 335
    end if

    ! case 2 - thin domain along Y
    if (ily==1) then
        nix(1:14)=nix4
        niy(1:14)=niy4
        niz(1:14)=niz4
        nsbcll=14
        if (ix==1) then
            nix(1:11)=nix3
            niy(1:11)=niy3
            niz(1:11)=niz3
            nsbcll=11
        elseif (ix==ilx) then
            nix(1:17)=nix2
            niy(1:17)=niy2
            niz(1:17)=niz2
            nsbcll=17
        end if    
        goto 335
    end if
    
    ! case 3 - normal domain   
    if (iy==1.and.ix<ilx) then
        nix(1:11)=nix3
        niy(1:11)=niy3
        niz(1:11)=niz3
        nsbcll=11
    elseif (iy==ily.and.ix>1) then
        nix(1:17)=nix2
        niy(1:17)=niy2
        niz(1:17)=niz2
        nsbcll=17
    else
        nix(1:14)=nix1
        niy(1:14)=niy1
        niz(1:14)=niz1
        nsbcll=14
    end if
335 continue
    
    ! put all data to arrays
    do kc=1,nsbcll
        jx=ix+nix(kc)
        jy=iy+niy(kc)
        jz=iz+niz(kc)
        cz=0
        if(jz>ilz)then
            jz=jz-ilz
            cz=dlz
        elseif(jz<1)then
            jz=jz+ilz
            cz=-dlz
        end if
        fmap(1,kc,ic)=jx
        fmap(2,kc,ic)=jy
        fmap(3,kc,ic)=jz
        fmapc(kc,ic)=cz
    end do
    
    ! move to the next subcell
    ix=ix+1
    if (ix>ilx) then
        ix=1
        iy=iy+1
        if(iy>ily)then
            iy=1
            iz=iz+1
        end if
    end if
end do

end

subroutine fcdrij (i,j,rsq,xd,yd,zd,arrtype)
! #####################################################################
! #                                                                   #
! #   subroutine 3:                                                   #
! #   calculate all forces: conservative, dissipative, and random     #
! #   for 1 pair of particles                                         #
! #                                                                   #
! #####################################################################
implicit none
save

common /datatech/ alpha,sigma,gamma
common /time/ dt
common /sizes/ sz
common /fxyz/ fx,fy,fz
common /vwxyz/ vwx,vwy,vwz
common /ptype/ kindp
common /commf/ fsend,rcou
common /prs/ press
common /reinit/ rein

real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) 
real*4,pointer :: fx(:),fy(:),fz(:)
integer*1,pointer :: kindp(:,:)
real*4,pointer :: fsend(:)
real*4,pointer :: alpha(:,:)

integer*4 sz
real*4 dt2inv,sigma,rfac,dt,gamma
real*4 rij,rsq,rijinv,exij,eyij,ezij,xd,yd,zd
real*4 omega
integer*4 rcou(4)
real*4 vxij,vyij,vzij
real*4 fcfac,fdfac,frfac,fctot,fcxij,fcyij,fczij
real*4 uni
integer*4 arrtype,nfi,nfj
integer*4 i,j
real*4 press(3)
logical init/.true./
integer*4 rein

! set initial values
if (init .or. rein==1) then
    init = .false.
    rfac=sqrt(3.d0) 
    dt2inv=sigma*rfac/sqrt(dt)
end if


! calculate distance between particles
rij=sqrt(rsq)  
rijinv=1.0/rij
omega=1.0 - rij

! calculate distance vector projections
exij=xd*rijinv
eyij=yd*rijinv
ezij=zd*rijinv

! calculate relative velocity
vxij=vwx(i,1) - vwx(j,arrtype)     
vyij=vwy(i,1) - vwy(j,arrtype)
vzij=vwz(i,1) - vwz(j,arrtype)

! get type of particles
nfi=kindp(i,1)             
nfj=kindp(j,arrtype)

! calculate the pairwise concervative force
fcfac=omega*alpha(nfi,nfj)

! calculate the pairwise dissipative force
fdfac=-gamma*omega*omega*(exij*vxij + eyij*vyij + ezij*vzij)
             
! calculate the pairwise random force
frfac=omega*dt2inv*(2.0*uni()-1.0)

! calculate the total force      
fctot =fcfac + fdfac + frfac

! calculate force projections
fcxij=fctot * exij
fcyij=fctot * eyij
fczij=fctot * ezij
fx(i)=fx(i) + fcxij
fy(i)=fy(i) + fcyij
fz(i)=fz(i) + fczij

! calculate pressure
press(3)=press(3)+fcfac*ezij*ezij/rijinv
press(2)=press(2)+fcfac*eyij*eyij/rijinv
press(1)=press(1)+fcfac*exij*exij/rijinv

! apply Newton 3rd law       
if (arrtype==1) then
    fx(j)=fx(j) - fcxij
    fy(j)=fy(j) - fcyij
    fz(j)=fz(j) - fczij
else
    ! apply forces to "ghost" beads
    fsend(3*j-2)=-fcxij+fsend(3*j-2)
    fsend(3*j-1)=-fcyij+fsend(3*j-1)
    fsend(3*j)=-fczij+fsend(3*j)
end if

end
      
subroutine fsendunpck()
! #####################################################################
! #                                                                   #
! #   subroutine 4:                                                   #
! #   send and receive forces                                         #
! #   send to the "right", receive from the "left"                    #
! #   update forces matrix acting on "real" particles                 #
! #                                                                   #
! #####################################################################
implicit none
include 'mpif.h'
save

common /funnumb/ npl,fann
common /fxyz/ fx,fy,fz
common /commf/ fsend,rcou
common /na/ natms,natmsfan
common /sndrcv/ sndr,sndl,fanl,fanr,rcv,rcvf
common /mpidpd/ rank, nproc
common /top/ cln,nrow,ncol

real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:) !size ssz,4;ssz,4;fsz,4;fsz,4;ssz*8,fsz*8
real*4,pointer :: fx(:),fy(:),fz(:) ! size sz 
real*4,pointer :: fsend(:) ! size 3*fsz/2
integer*4,pointer :: fann(:,:) !size fsz/8,4
integer*4 npl(4),natms,natmsfan,i,j,k,m,indcur
integer*4 rcou(4)
integer*4 rank, nproc,cln,nrow,ncol
integer lproc(4),rproc(4)
integer ierr
integer status(MPI_STATUS_SIZE)
logical init /.true./

! set nodes exchange matrix
if (init) then
    lproc(3)=rank-1
    if (ncol==1) lproc(3)=rank+cln-1
    lproc(2)=lproc(3)-cln
    if (nrow==1) lproc(2)=nproc-cln+lproc(3)
    lproc(4)=lproc(3)+cln
    if (nrow==nproc/cln) lproc(4)=lproc(3)-nproc+cln
    lproc(1)=rank-cln
    if (nrow==1) lproc(1)=nproc-cln+rank
    rproc(3)=rank+1
    if (ncol==cln) rproc(3)=rank-cln+1
    rproc(2)=rproc(3)+cln
    if (nrow==nproc/cln) rproc(2)=rproc(3)-nproc+cln
    rproc(4)=rproc(3)-cln
    if (nrow==1) rproc(4)=nproc-cln+rproc(3)
    rproc(1)=rank+cln
    if (nrow==nproc/cln) rproc(1)=rank-nproc+cln
    init=.false.
end if

! send and receive forces; update forces matrix
indcur=1
do m=1,4
    if (rproc(m)==rank) then
        rcvf(1:rcou(m))=fsend(indcur:indcur-1+rcou(m))
    else
        call MPI_Sendrecv(fsend(indcur), rcou(m), MPI_real,rproc(m), 1,rcvf, 3*npl(m), MPI_real, lproc(m), 1, MPI_Comm_world, status,ierr )
    end if
    indcur=indcur+rcou(m)
    j=1
    do i=1,3*npl(m),3
        k=fann(j,m)
        fx(k)=fx(k)+rcvf(i)
        fy(k)=fy(k)+rcvf(i+1)
        fz(k)=fz(k)+rcvf(i+2)
        j=j+1
    end do
end do
fsend=0

end

real*4 function uni()
! #####################################################################
! #                                                                   #
! #   function 5:                                                     #
! # random number generator based on the universal random number      #
! # generator of marsaglia, zaman and tsang                           #
! # (stats and prob. lett. 8 (1990) 35-39.)                           #
! #                                                                   #
! # it must be called once to initialise parameters u,c,cd,cm         #
! #                                                                   #
! #####################################################################
implicit none

common /mpidpd/ rank, nproc
logical,save :: init
data init /.true./
integer,save :: ir,jr
integer i,ii,j,jj,k,l,m,idnode
real*8,save :: c,cd,cm,u(1:97)
real*8 s,t
integer*4 nproc,rank
      

! initialise parameters u,c,cd,cm
if (init) then
    init = .false.

! initial values of i,j,k must be in range 1 to 178 (not all 1)
! initial value of l must be in range 0 to 168
    i = mod(rank,166) + 12
    j = mod(rank,144) + 34
    k = mod(rank,122) + 56
    l = mod(rank,90)  + 78
    ir = 97
    jr = 33
    do ii=1,97
        s = 0.0d0
        t = 0.5d0
        do jj=1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if (mod(l*m,64)>=32) s = s+t
            t = 0.5d0*t
        end do
        u(ii)=s
    end do
    c  =   362436.0d0/16777216.0d0
    cd =  7654321.0d0/16777216.0d0
    cm = 16777213.0d0/16777216.0d0
end if

! calculate random number
uni=u(ir)-u(jr)
if (uni < 0.0d0) uni = uni + 1.0d0
u(ir)=uni
ir=ir-1
if (ir == 0) ir = 97
jr=jr-1
if (jr == 0) jr = 97
c = c-cd
if (c < 0.0d0) c = c+cm
uni = uni-c
if (uni < 0.0d0) uni = uni + 1.0d0

end

 subroutine intr()
! #################################################################
! #                                                               #
! #    subroutine 6:                                              #
! #    integrate the positions of the particles,                  #
! #    send leaving and ghost particles                           #
! #                                                               #
! #################################################################
implicit none
save

common /sizes/ sz
common /time/ dt
common /na/ natms,natmsfan
common /nab/ natmsfanbond
common /comm/ ssz,fsz
common /npall/ npall
common /rxyz/ rx,ry,rz
common /fxyz/ fx,fy,fz
common /vxyz/ vx,vy,vz
common /vwxyz/ vwx,vwy,vwz
common /correl/ corr,bcorr,bcorrt
common /ptype/ kindp
common /cell/ dlx,dly,dlz,rcut,kpbc
common /sndrcv/ sndr,sndl,fanl,fanr,rcv,rcvf
common /funnumb/ npl,fann
common /shwi/ shwi
common /lambd/ lambda
common /reinit/ rein

integer*4 sz,ssz,fsz,npall

real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:)
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 
real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) 
real*4,pointer :: fx(:),fy(:),fz(:)
real*4,pointer :: vx(:),vy(:),vz(:) 
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:)
integer*1,pointer :: kindp(:,:) 
integer*4,pointer :: fann(:,:)

real*4 dt
real*4 vwxb,vwyb,vwzb,vxb,vyb,vzb,rxb,ryb,rzb
real*4 lambda,dtlam,dt2
real*4 dlx,dly,dlz,rcut,shwi
integer*4 natms,natmsfan,natmsfanbond
integer*4 i,j,shft,r1(4),l1(4),lr(4),rl(4)
integer*4 res,fan,kpbc,resfan,neinum
integer*4 npl(4)
integer*4 neinumb(2),neilist(4,2)
real*4 rxres(4,2),ryres(4,2)
logical init/.true./
integer*4 rein

! set initial values
if (init .or. rein==1) then
      dtlam = lambda*dt
      dt2 = 0.50*dt
      init=.false.
end if

! zero arrays
shft=0
r1=0
l1=0
lr=0
rl=0
fann=0
natmsfan=0
natmsfanbond=0
npl=0
bcorrt=0

! loop over particles, integrate positions and velocities 
do i = 1, natms
    vwxb = vx(i-shft) + dtlam*fx(i-shft)
    vwyb = vy(i-shft) + dtlam*fy(i-shft)
    vwzb = vz(i-shft) + dtlam*fz(i-shft)
    vxb  = vx(i-shft) + dt2*fx(i-shft)
    vyb  = vy(i-shft) + dt2*fy(i-shft)
    vzb  = vz(i-shft) + dt2*fz(i-shft)
    if (dt*vxb>dlx) call error(60)
    if (dt*vyb>dly) call error(61)
    rxb  = rx(i-shft,1) + dt*vxb
    ryb  = ry(i-shft,1) + dt*vyb
    rzb  = rz(i-shft,1) + dt*vzb

    ! check if ith bead will leave this node's subvolume after integration
    call pbc(rxb,ryb,rzb,res,neinum)

    ! put information about leaving bead to the communication array 
    if (res==1) then
        sndr(r1(neinum)+1,neinum)=rxb
        sndr(r1(neinum)+2,neinum)=ryb
        sndr(r1(neinum)+3,neinum)=rzb
        sndr(r1(neinum)+4,neinum)=vxb
        sndr(r1(neinum)+5,neinum)=vyb
        sndr(r1(neinum)+6,neinum)=vzb
        sndr(r1(neinum)+7,neinum)=vwxb
        sndr(r1(neinum)+8,neinum)=vwyb
        sndr(r1(neinum)+9,neinum)=vwzb
        sndr(r1(neinum)+10,neinum)=real(corr(i-shft,1))
        r1(neinum)=r1(neinum)+10
        vx(i-shft)=vx(natms-shft)
        vy(i-shft)=vy(natms-shft)
        vz(i-shft)=vz(natms-shft)
        fx(i-shft)=fx(natms-shft)
        fy(i-shft)=fy(natms-shft)
        fz(i-shft)=fz(natms-shft)
        rx(i-shft,1)=rx(natms-shft,1)
        ry(i-shft,1)=ry(natms-shft,1)
        rz(i-shft,1)=rz(natms-shft,1)
        kindp(i-shft,1)=kindp(natms-shft,1)
        corr(i-shft,1)=corr(natms-shft,1)
        shft=shft+1
    elseif (res==-1) then
        sndl(l1(neinum)+1,neinum)=rxb
        sndl(l1(neinum)+2,neinum)=ryb
        sndl(l1(neinum)+3,neinum)=rzb
        sndl(l1(neinum)+4,neinum)=vxb
        sndl(l1(neinum)+5,neinum)=vyb
        sndl(l1(neinum)+6,neinum)=vzb
        sndl(l1(neinum)+7,neinum)=vwxb
        sndl(l1(neinum)+8,neinum)=vwyb
        sndl(l1(neinum)+9,neinum)=vwzb
        sndl(l1(neinum)+10,neinum)=real(corr(i-shft,1))
        l1(neinum)=l1(neinum)+10
        vx(i-shft)=vx(natms-shft)
        vy(i-shft)=vy(natms-shft)
        vz(i-shft)=vz(natms-shft)
        fx(i-shft)=fx(natms-shft)
        fy(i-shft)=fy(natms-shft)
        fz(i-shft)=fz(natms-shft)
        rx(i-shft,1)=rx(natms-shft,1)
        ry(i-shft,1)=ry(natms-shft,1)
        rz(i-shft,1)=rz(natms-shft,1)
        kindp(i-shft,1)=kindp(natms-shft,1)
        corr(i-shft,1)=corr(natms-shft,1)
        shft=shft+1

    ! update positions and velocities arrays if the bead remains inside node's subvolume
    else
        vx(i-shft)=vxb
        vy(i-shft)=vyb
        vz(i-shft)=vzb
        rx(i-shft,1)=rxb
        ry(i-shft,1)=ryb
        rz(i-shft,1)=rzb
        vwx(i-shft,1)=vwxb
        vwy(i-shft,1)=vwyb
        vwz(i-shft,1)=vwzb
    end if
end do

! calculate remaining number of beads; zero arrays
natms=natms-shft
vx(natms+1:)=0
vy(natms+1:)=0
vz(natms+1:)=0
rx(natms+1:,1)=0
ry(natms+1:,1)=0
rz(natms+1:,1)=0
fx(natms+1:)=0
fy(natms+1:)=0
fz(natms+1:)=0
vwx(natms+1:,1)=0
vwy(natms+1:,1)=0
vwz(natms+1:,1)=0

! send and receive leaving atoms
call send(l1,r1,1)

! update beads arrays
call unpck(l1,r1,1)

! define ghost beads
do i=1,natms
    rxb=rx(i,1)
    ryb=ry(i,1)
    call pbcf(rxb,ryb,neinumb,neilist,rxres,ryres)
    do j=1,neinumb(1)
        neinum=neilist(j,1)
        npl(neinum)=npl(neinum)+1
        fann(npl(neinum),neinum)=i
        fanl(lr(neinum)+1,neinum)=rxres(j,1)
        fanl(lr(neinum)+2,neinum)=ryres(j,1)
        fanl(lr(neinum)+3,neinum)=rz(i,1)
        fanl(lr(neinum)+4,neinum)=vwx(i,1)
        fanl(lr(neinum)+5,neinum)=vwy(i,1)
        fanl(lr(neinum)+6,neinum)=vwz(i,1)
        fanl(lr(neinum)+7,neinum)=real(corr(i,1))
        lr(neinum)=lr(neinum)+7
    end do
    
    do j=1,neinumb(2)
        neinum=neilist(j,2)
        fanr(rl(neinum)+1,neinum)=rxres(j,2)
        fanr(rl(neinum)+2,neinum)=ryres(j,2)
        fanr(rl(neinum)+3,neinum)=rz(i,1)
        fanr(rl(neinum)+4,neinum)=real(corr(i,1))
        rl(neinum)=rl(neinum)+4
    end do
end do

! send and receive ghost beads
call send(lr,rl,2)

! update ghost beads arrays
call unpck(lr,rl,2)

! correct arrays of global bead numbers 
do i=1,natms
    bcorrt(corr(i,1))=1
    bcorr(corr(i,1))=i
end do

end
      
subroutine pbc(rxb,ryb,rzb,res,neinum)
! #####################################################################
! #                                                                   #
! #   subroutine 7:                                                   #
! #   apply periodic boundary conditions                              #
! #   define destination for leaving beads                            #
! #                                                                   #
! #####################################################################
implicit none
     
common /cell/ dlx,dly,dlz,rcut,kpbc
common /shwi/ shwi

integer*4 kpbc
real*4 dlx,dly,dlz,rcut
real*4 shwi
integer*4 res,neinum
real*4 rxb,ryb,rzb
rzb=rzb-dlz*nint(rzb/dlz-0.5)
res=0
neinum=0

! check which boundary was crossed
if (rxb>dlx) then
    rxb=rxb-dlx
    res=1
    if (ryb>dly) then
        neinum=2
        ryb=ryb-dly
    elseif (ryb<0) then
        neinum=4
        ryb=ryb+dly
    else 
        neinum=3
    end if
    return       
elseif (rxb<0) then
    rxb=rxb+dlx
    res=-1
    if (ryb>dly) then
        neinum=4
        ryb=ryb-dly
    elseif (ryb<0) then
        neinum=2
        ryb=ryb+dly
    else 
        neinum=3
    end if
    return
end if
if (ryb>dly) then
    ryb=ryb-dly
    res=1
    neinum=1
elseif (ryb<0) then
    ryb=ryb+dly
    res=-1
    neinum=1
end if

end

subroutine pbcf(rxb,ryb,neinumb,neilist,rxres,ryres)
! #####################################################################
! #                                                                   #
! #   subroutine 8:                                                   #
! #   apply PBC                                                       #
! #   define destination for ghost atoms                              #
! #                                                                   #
! #####################################################################
implicit none

common /cell/ dlx,dly,dlz,rcut,kpbc
common /shwi/ shwi

integer*4 kpbc
real*4 dlx,dly,dlz,rcut
real*4 shwi
real*4 rxb,ryb
integer*4 neinumb(2),neilist(4,2)
real*4 rxres(4,2),ryres(4,2)

! set initial values
neinumb=0
neilist=0
rxres=rxb
ryres=ryb

! check which node should receive information about current bead
if (rxb<shwi) then
    neinumb(1)=neinumb(1)+1
    neilist(neinumb(1),1)=3
    rxres(neinumb(1),1)=rxb+dlx
    
    if(ryb<shwi) then
        neinumb(1)=neinumb(1)+1
        neilist(neinumb(1),1)=2
        ryres(neinumb(1),1)=ryb+dly
        rxres(neinumb(1),1)=rxb+dlx
    end if
    if(ryb>dly-shwi) then
        neinumb(1)=neinumb(1)+1
        neilist(neinumb(1),1)=4
        ryres(neinumb(1),1)=ryb-dly
        rxres(neinumb(1),1)=rxb+dlx
    end if
end if

if (rxb>dlx-shwi) then
    neinumb(2)=neinumb(2)+1
    neilist(neinumb(2),2)=3
    rxres(neinumb(2),2)=rxb-dlx
    if(ryb<shwi) then
        neinumb(2)=neinumb(2)+1
        neilist(neinumb(2),2)=4
        ryres(neinumb(2),2)=ryb+dly
        rxres(neinumb(2),2)=rxb-dlx
    end if
    if(ryb>dly-shwi) then
        neinumb(2)=neinumb(2)+1
        neilist(neinumb(2),2)=2
        ryres(neinumb(2),2)=ryb-dly
        rxres(neinumb(2),2)=rxb-dlx
    end if
end if

if(ryb>dly-shwi) then
    neinumb(2)=neinumb(2)+1
    neilist(neinumb(2),2)=1
    ryres(neinumb(2),2)=ryb-dly
        
end if
if(ryb<shwi) then
    neinumb(1)=neinumb(1)+1
    neilist(neinumb(1),1)=1
    ryres(neinumb(1),1)=ryb+dly
        
end if

end

subroutine send(l1,r1,flag)
! #####################################################################
! #                                                                   #
! #   subroutine 9:                                                   #
! #   send and receive information about beads                        #
! #   flag=1 for leaving beads                                        #
! #   flag=2 for ghost beads                                          #
! #                                                                   #
! #####################################################################
implicit none
save
include "mpif.h"

common /comm/ ssz,fsz
common /mpidpd/ rank, nproc
common /top/ cln,nrow,ncol 
common /sndrcv/ sndr,sndl,fanl,fanr,rcv,rcvf

real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:) 

integer*4 ssz,fsz
integer*4 rank, nproc,cln,nrow,ncol
integer*4 r1(4),l1(4)

integer*4 r(4),l(4)
integer*4 i,indcur
integer*1 flag
integer lproc(4),rproc(4)
integer ierr
integer status(MPI_STATUS_SIZE)
logical init /.true./

! construct nodes interaction matrix
if (init) then
    lproc(3)=rank-1
    if (ncol==1) lproc(3)=rank+cln-1
    lproc(2)=lproc(3)-cln
    if (nrow==1) lproc(2)=nproc-cln+lproc(3)
    lproc(4)=lproc(3)+cln
    if (nrow==nproc/cln) lproc(4)=lproc(3)-nproc+cln
    lproc(1)=rank-cln
    if (nrow==1) lproc(1)=nproc-cln+rank
    rproc(3)=rank+1
    if (ncol==cln) rproc(3)=rank-cln+1
    rproc(2)=rproc(3)+cln
    if (nrow==nproc/cln) rproc(2)=rproc(3)-nproc+cln
    rproc(4)=rproc(3)-cln
    if (nrow==1) rproc(4)=nproc-cln+rproc(3)
    rproc(1)=rank+cln
    if (nrow==nproc/cln) rproc(1)=rank-nproc+cln
    init=.false.
end if

! send and receive number of beads
do i=1,4
    if (rproc(i)==rank) then
        l(i)=r1(i)
    else
        call MPI_Sendrecv( r1(i), 1, MPI_integer,rproc(i), 1,l(i), 1, MPI_integer, lproc(i), 1, MPI_Comm_world, status,ierr )
    end if
end do

do i=1,4
    if (rproc(i)==rank) then
        r(i)=l1(i)
    else
        call MPI_Sendrecv( l1(i), 1, MPI_integer,lproc(i), 1,r(i), 1, MPI_integer, rproc(i), 1, MPI_Comm_world, status,ierr )
    end if
end do

! send and receive information about leaving beads
if (flag==1) then
    indcur=1
    do i=1,4
        if (rproc(i)==rank) then
            rcv(indcur:indcur+l(i)-1)=sndr(1:l(i),i)
        else
            call MPI_Sendrecv( sndr(:,i),r1(i),MPI_real,rproc(i),1,rcv(indcur),l(i),MPI_real,lproc(i),1,MPI_Comm_world,status,ierr )
        end if
        indcur=indcur+l(i)
    end do
    do i=1,4
        if (rproc(i)==rank) then
            rcv(indcur:indcur+r(i)-1)=sndl(1:r(i),i)
        else
            call MPI_Sendrecv( sndl(:,i),l1(i),MPI_real,lproc(i),1,rcv(indcur),r(i),MPI_real,rproc(i),1,MPI_Comm_world,status,ierr )
        end if
        indcur=indcur+r(i)
    end do
    r1=r
    l1=l
end if

! send and receive information about ghost beads
if (flag==2) then
    indcur=1
    do i=1,4
        if (rproc(i)==rank) then
            rcvf(indcur:indcur+r(i)-1)=fanl(1:r(i),i)
        else
            call MPI_Sendrecv( fanl(:,i),l1(i),MPI_real,lproc(i),1,rcvf(indcur),r(i),MPI_real,rproc(i),1,MPI_Comm_world,status,ierr )
        end if
        indcur=indcur+r(i)
    end do
    do i=1,4
        if (rproc(i)==rank) then
            rcvf(indcur:indcur+l(i)-1)=fanr(1:l(i),i)
        else
            call MPI_Sendrecv( fanr(:,i),r1(i),MPI_real,rproc(i),1,rcvf(indcur),l(i),MPI_real,lproc(i),1,MPI_Comm_world,status,ierr )
        end if
        indcur=indcur+l(i)
    end do
    r1=r
    l1=l
end if

end

subroutine unpck(l1,r1,flag)
! #####################################################################
! #                                                                   #
! #   subroutine 10:                                                  #
! #   unpack data, update arrays with new beads information           #
! #   flag=1 for leaving atoms                                        #
! #   flag=2 for ghost atoms                                          #
! #                                                                   #
! #####################################################################
implicit none

common /na/ natms,natmsfan
common /nab/ natmsfanbond
common /comm/ ssz,fsz
common /npall/ npall
common /sizes/ sz
common /rxyz/ rx,ry,rz
common /vxyz/ vx,vy,vz
common /vwxyz/ vwx,vwy,vwz
common /correl/ corr,bcorr,bcorrt
common /ptype/ kindp
common /sndrcv/ sndr,sndl,fanl,fanr,rcv,rcvf
common /ptypeall/ kindpall !particle types
common /cell/ dlx,dly,dlz,rcut,kpbc
common /shwi/ shwi
common /commf/ fsend,rcou

real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:) !size ssz,4;ssz,4;fsz,4;fsz,4;ssz*8,fsz*8

integer*4 natms,natmsfan,natmsfanbond
integer*4 sz,ssz,fsz,npall
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) ! size sz,3
real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) ! size sz,2
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
integer*4,pointer :: corr(:,:),bcorr(:)! size sz,3 and npall
integer*1,pointer :: bcorrt(:) !size npall
integer*1,pointer :: kindp(:,:) ! size sz,3
integer*1,pointer :: kindpall(:) ! size npall
real*4,pointer :: fsend(:) ! size 3*fsz/2

real*4 dlx,dly,dlz,rcut,shwi
integer*4 r1(4),l1(4),rcou(4)

integer*4 rcvb
real*4 rxb
integer*4 i
integer*1 kpbc,flag

if (flag==1) then

    ! update arrays of real beads
    do i=0,sum(l1)+sum(r1)-1,10
        rx(natms+1,1)=rcv(i+1)   
        ry(natms+1,1)=rcv(i+2)
        rz(natms+1,1)=rcv(i+3)
        vx(natms+1)=rcv(i+4)
        vy(natms+1)=rcv(i+5)
        vz(natms+1)=rcv(i+6)
        vwx(natms+1,1)=rcv(i+7)
        vwy(natms+1,1)=rcv(i+8)
        vwz(natms+1,1)=rcv(i+9)
        rcvb= int(rcv(i+10))
        corr(natms+1,1)=rcvb
        bcorr(rcvb)=natms+1
        bcorrt(rcvb)=1
        kindp(natms+1,1)=kindpall(rcvb)
        natms=natms+1
    end do
else

    ! update arrays of ghost beads
    rcou=r1/7*3
    do i=0,sum(r1)-1,7
        rx(natmsfan+1,2)=rcvf(i+1)
        ry(natmsfan+1,2)=rcvf(i+2)
        rz(natmsfan+1,2)=rcvf(i+3)
        vwx(natmsfan+1,2)=rcvf(i+4)
        vwy(natmsfan+1,2)=rcvf(i+5)
        vwz(natmsfan+1,2)=rcvf(i+6)
        rcvb=int(rcvf(i+7))
        corr(natmsfan+1,2)=rcvb
        bcorr(rcvb)=natmsfan+1
        bcorrt(rcvb)=2
        kindp(natmsfan+1,2)=kindpall(rcvb)
        natmsfan=natmsfan+1
    end do
    do i=sum(r1),sum(l1)+sum(r1)-1,4
        rx(natmsfanbond+1,3)=rcvf(i+1)
        ry(natmsfanbond+1,3)=rcvf(i+2)
        rz(natmsfanbond+1,3)=rcvf(i+3)
        rcvb=int(rcvf(i+4))
        corr(natmsfanbond+1,3)=rcvb
        bcorr(rcvb)=natmsfanbond+1
        bcorrt(rcvb)=3
        kindp(natmsfanbond+1,3)=kindpall(rcvb)
        natmsfanbond=natmsfanbond+1
    end do
end if

end   

subroutine initial()
! #################################################################
! #                                                               #
! #    subroutine 11:                                             #
! #    create initial ghost beads arrays                          #
! #                                                               #
! #################################################################
implicit none

common /rxyz/ rx,ry,rz
common /vxyz/ vx,vy,vz
common /vwxyz/ vwx,vwy,vwz
common /fxyz/ fx,fy,fz
common /ptype/ kindp
common /correl/ corr,bcorr,bcorrt
common /den/ rho
common /sizes/ sz
common /npall/ npall
common /na/ natms,natmsfan
common /nab/ natmsfanbond
common /cell/ dlx,dly,dlz,rcut,kpbc
common /mpidpd/ rank, nproc
common /shwi/ shwi
common /comm/ ssz,fsz
common /sndrcv/ sndr,sndl,fanl,fanr,rcv,rcvf
common /funnumb/ npl,fann

real*4 uni
integer*4 fsz,ssz
integer*4 sz,npall
real*4 dlx,dly,dlz,rcut,rho,vall(3),shwi,vl
integer*4 neinumb(2),neilist(4,2)
real*4 rxres(4,2),ryres(4,2),rxb,ryb
integer*4 kpbc,fan,res,resfan
integer*4 natms,natmsfan,natmsfanbond
integer*4 i,lr(4),npl(4),rl(4),j,neinum
integer*4 rank,nproc
real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:) !size ssz,4;ssz,4;fsz,4;fsz,4;ssz*8,fsz*8
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) ! size sz,3
real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) ! size sz,2
real*4,pointer :: fx(:),fy(:),fz(:) ! size sz
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
integer*1,pointer :: kindp(:,:) ! size sz,3
integer*4,pointer :: corr(:,:),bcorr(:)! size sz,3 and npall
integer*1,pointer :: bcorrt(:) !size npall
integer*4,pointer :: fann(:,:) !size fsz/8,4

! set initial values
vwx(:,1)=vx
vwy(:,1)=vy
vwz(:,1)=vz
fx=0
fy=0
fz=0
lr=0
rl=0
fann=0
npl=0

! define ghost beads
do i=1,natms
    rxb=rx(i,1)
    ryb=ry(i,1)
    call pbcf(rxb,ryb,neinumb,neilist,rxres,ryres)
    do j=1,neinumb(1)
        neinum=neilist(j,1)
        npl(neinum)=npl(neinum)+1
        fann(npl(neinum),neinum)=i
        fanl(lr(neinum)+1,neinum)=rxres(j,1)
        fanl(lr(neinum)+2,neinum)=ryres(j,1)
        fanl(lr(neinum)+3,neinum)=rz(i,1)
        fanl(lr(neinum)+4,neinum)=vwx(i,1)
        fanl(lr(neinum)+5,neinum)=vwy(i,1)
        fanl(lr(neinum)+6,neinum)=vwz(i,1)
        fanl(lr(neinum)+7,neinum)=real(corr(i,1))
        lr(neinum)=lr(neinum)+7
    end do
    
    do j=1,neinumb(2)
        neinum=neilist(j,2)
        fanr(rl(neinum)+1,neinum)=rxres(j,2)
        fanr(rl(neinum)+2,neinum)=ryres(j,2)
        fanr(rl(neinum)+3,neinum)=rz(i,1)
        fanr(rl(neinum)+4,neinum)=real(corr(i,1))
        rl(neinum)=rl(neinum)+4
    end do
end do
natmsfan=0
natmsfanbond=0

! send and receive ghost beads
call send(lr,rl,2)

! update ghost beads arrays
call unpck(lr,rl,2)

! correct arrays of global bead numbers 
do i=1,natms
    bcorrt(corr(i,1))=1
    bcorr(corr(i,1))=i
end do

end

subroutine intv()     
! #################################################################
! #                                                               #
! #    subroutine 13:                                             #
! #    integrate velocities using a time step of "h"              #
! #                                                               #
! #################################################################
implicit none
save

common /fxyz/ fx,fy,fz
common /vxyz/ vx,vy,vz
common /time/ dt
common /na/ natms,natmsfan
common /reinit/ rein
real*4,pointer :: fx(:),fy(:),fz(:) ! size sz
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
real*4 dt,dt2
integer*4 i,natms,natmsfan
logical init/.true./
integer*4 rein

! set initial value
if (init .or. rein==1) then
    dt2 = 0.50*dt
    init=.false.
    rein=0
end if
      
! integrate velocities
do i=1, natms
    vx(i) = vx(i) + dt2*fx(i)
    vy(i) = vy(i) + dt2*fy(i)
    vz(i) = vz(i) + dt2*fz(i)
end do

end

subroutine readrst()
! #################################################################
! #                                                               #
! #    subroutine 15:                                             #
! #    read restart file                                          #
! #                                                               #
! #################################################################
implicit none

common /cell/ dlx,dly,dlz,rcut,kpbc
common /cella/ dlxa,dlya,dlza
common /npall/ npall
common /mpidpd/ rank, nproc
common /na/ natms,natmsfan
common /rxyz/ rx,ry,rz
common /vxyz/ vx,vy,vz
common /bxsize/ ilx,ily,ilz
common /sizes/ sz
common /den/ rho 
common /comm/ ssz,fsz
common /shwi/ shwi
common /ptype/ kindp
common /correl/ corr,bcorr,bcorrt
common /bnd/ cn
common /rst/ nbond
common /typel/ typelist,ntype
common /ptypeall/ kindpall
common /top/ cln,nrow,ncol 
common /vale/ val

integer*4 natms,natmsfan,count
real*4 dlx,dly,dlz,rcut
real*4 dlxa,dlya,dlza,cell(9)
integer*4 npall,nb,nc,kpbc
integer*4 rank,nproc,cln,nrow,ncol
integer*4 ilx,ily,ilz
real*4 rho
real*4 shwi

integer*4 sz,ssz,fsz
integer*4 corrb
integer*1 kindpb,valb
real*4 rxb,ryb,rzb,vxb,vyb,vzb
integer*4 i,j,k,i1,i2,i3,bt,at
logical exists
character*9   dummy
integer*4 nbond,kb
integer*1 typelist(10),ntype

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 
real*4,pointer :: vx(:),vy(:),vz(:) 
integer*1,pointer :: kindp(:,:) !
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:) 
integer*4,pointer :: cn(:,:) 
integer*1,pointer :: kindpall(:) 
integer*4,pointer :: bond(:) 
integer*1,pointer :: val(:) 

inquire( file = 'restart.dat', exist = exists )
if ( .not. exists ) call error(153)
open(1, file = 'restart.dat', status = 'old' )

! read number of beads and box sizes
read(1,*, err=10, end=10) npall,rho
read(1,*, err=10, end=10) dlxa,dlya,dlza

! calculate nodes topology
rcut=1.0
cln=int(sqrt(real(nproc)))
do
    if (mod(nproc,cln)==0) exit
    cln=cln+1
end do
nrow=rank/cln+1
ncol=mod(rank,cln)+1

! calculate box size on each node
dlz=dlza
dly=real(dlya/real(nproc/cln))
if (dly<rcut) call error(150)
dlx=real(dlxa/real(cln))
if (dlx<rcut) call error(151)

! calculate number of subcells
ilx=int(dlx/rcut)
ily=int(dly/rcut)
ilz=int(dlz/rcut)

! sizes of arrays
sz=2*max0(int(dlz*dly*dlx*rho*3),int(3*(2*dlz*rho*dlx*shwi+2*dlz*rho*dly*shwi+4*dlz*rho*shwi*shwi)))
ssz=2*max0(int(shwi*dly*dlz*11*rho),int(shwi*dlx*dlz*11*rho))
fsz=2*max0(int(shwi*dly*dlz*8*3*rho),int(shwi*dlx*dlz*8*3*rho))

! allocate arrays
allocate (rx(sz,3),ry(sz,3),rz(sz,3),vx(sz),vy(sz),vz(sz))
allocate (kindp(sz,3))
allocate (corr(sz,3),bcorr(npall),bcorrt(npall))
allocate (cn(npall,0:0))
allocate (kindpall(npall))
allocate (val(npall))

! zero arrays
natms=0
corr=0
bcorr=0
bcorrt=0
rx=0
ry=0
rz=0
vx=0
vy=0
vz=0
cn=0
ntype=0
typelist=0

! read information about beads
! add fix: type 1 - chain, boxsize = dlxa/2-1; type 2 - solvent, boxsize = dlxa
do i = 1, npall
    read(1,'(i8,2i4,3f14.6)', err=10, end=10)corrb,valb,kindpb,rxb,ryb,rzb

        ! move current bead inside the simulation box
    !if (kindpb==1) then
     !   if (rxb>(dlxa/2-1)) rxb=rxb-int(rxb/(dlxa/2-1))*(dlxa/2-1)
      !  if (rxb<0) rxb=rxb+(dlxa/2-1)*(abs(int(rxb/(dlxa/2-1)))+1)
       ! if (ryb>(dlya/2-1)) ryb=ryb-int(ryb/(dlya/2-1))*(dlya/2-1)
        !if (ryb<0) ryb=ryb+(dlya/2-1)*(abs(int(ryb/(dlya/2-1)))+1)
        !if (rzb>(dlza/2-1)) rzb=rzb-int(rzb/(dlza/2-1))*(dlza/2-1)
        !if (rzb<0) rzb=rzb+(dlza/2-1)*(abs(int(rzb/(dlza/2-1)))+1)
        !kindpall(corrb)=kindpb
        !val(corrb)=valb
    if (kindpb==1 .or. kindpb==2) then
        if (rxb>dlxa) rxb=rxb-int(rxb/dlxa)*dlxa
        if (rxb<0) rxb=rxb+dlxa*(abs(int(rxb/dlxa))+1)
        if (ryb>dlya) ryb=ryb-int(ryb/dlya)*dlya
        if (ryb<0) ryb=ryb+dlya*(abs(int(ryb/dlya))+1)
        if (rzb>dlza) rzb=rzb-int(rzb/dlza)*dlza
        if (rzb<0) rzb=rzb+dlza*(abs(int(rzb/dlza))+1)
        kindpall(corrb)=kindpb
        val(corrb)=valb
    else
        stop 'Type of the bead does not equal 1 or 2. Check the restart.dat file.'
    end if
  
    ! obtain list of avaliable beads types
    do j = 1, ntype
        if ( kindpb == typelist(j) ) goto 1
    end do
    ntype = ntype + 1
    typelist(ntype) = kindpb
1   continue

    ! define node which current bead belongs to
    if ((rxb<ncol*dlx+0.00001*int(ncol/cln)).and.(rxb>=(ncol-1)*dlx).and.(ryb<nrow*dly+0.00001*int(nrow*cln/nproc)).and.(ryb>=(nrow-1)*dly)) then
        rx(natms+1,1)=rxb-(ncol-1)*dlx
        ry(natms+1,1)=ryb-(nrow-1)*dly
        rz(natms+1,1)=rzb
        kindp(natms+1,1)=kindpb
        corr(natms+1,1)=corrb
        bcorr(corrb)=natms+1
        bcorrt(corrb)=1
        natms=natms+1
    end if       
end do
rxb=natms ! rxb is just a buffer for summation
call gsum(rxb,1)

! check if all beads are inside the simulation box
if (rxb/=npall) call error(152)
count=0

! read bonds information
read(1,'(a9,i9)', err=10, end=10) dummy,nbond
allocate(bond(nbond*2))
kb = 1
do i = 1, nbond
    read(1,*, err=10, end=10) i1, i2
    bond(kb)=i1
    bond(kb+1)=i2
    kb=kb+2
    cn(i1,0) = cn(i1,0) + 1
    cn(i2,0) = cn(i2,0) + 1
end do

! define maximum number of bonds per bead
cn(:,0)=cn(:,0)+val(:)
kb=maxval(cn(:,0))
nullify (cn)

! allocate connection matrix array using maximum number of bonds per bead
allocate (cn(npall,0:kb))
cn=0
kb=1

! generate connection matrix array
do i = 1, nbond
    i1=bond(kb)
    i2=bond(kb+1)
    kb=kb+2
    cn(i1,0) = cn(i1,0) + 1
    cn(i1,cn(i1,0)) = i2
    cn(i2,0) = cn(i2,0) + 1
    cn(i2,cn(i2,0)) = i1
end do

! delete temporary bonds array
nullify(bond)

close(1)

! set initial velocities to 0
vx=0
vy=0
vz=0
return

10 call error(154)

end

subroutine forcesbond()
! #################################################################
! #                                                               #
! #    subroutine 16:                                             #
! #    calculate bond forces                                      #
! #                                                               #
! #################################################################
implicit none
include "mpif.h"
save

common /bondf/ k_eq,k_bond
common /rxyz/ rx,ry,rz
common /fxyz/ fx,fy,fz
common /na/ natms,natmsfan
common /cella/ dlxa,dlya,dlza
common /correl/ corr,bcorr,bcorrt
common /bnd/ cn
common /mpidpd/ rank, nproc
common /sizes/ sz
common /cell/ dlx,dly,dlz,rcut,kpbc
common /top/ cln,nrow,ncol 
common /prs/ press 
common /ptypeall/ kindpall
common /basndrcfa/ sndb,rcvb
common /bmult/ mult
common /rst/ nbond
common /shwi/ shwi

real*4 dlx,dly,dlz,rcut
integer*4 natms,natmsfan,mbnum,allnum
integer nam,ierr
real*4 dlxa,dlya,dlza
integer*4 i,i1,i2,i3,j,nb,k,i2b,shft,mult
integer*1 arrtype,flag
integer*4 rank, nproc,cln,nrow,ncol
real*4 rxij,ryij,rzij,rijsq
integer*4 sz
logical init/.true./
integer*4 kpbc,nbond
integer*4 arrsize
real*4 press(3)
real*4 shwi

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 
integer*4,pointer :: mb(:,:) 
real*4,pointer :: fx(:),fy(:),fz(:) 
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:) 
integer*4,pointer :: cn(:,:) 
integer,pointer :: rcount(:),displ(:) 
real*4,pointer :: sndb(:), rcvb(:) 
real*4,pointer :: k_eq(:,:),k_bond(:,:) 
integer*1,pointer :: kindpall(:) 

if (nbond==0) return

! allocate send/recieve arrays for bonds calculation
if (init) then
    init=.false.
    arrsize=mult*int(18*(2*dlx*dlz+2*dly*dlz)/(shwi**2))
    allocate( mb(sz,2),rcount(nproc),displ(nproc),sndb(arrsize),rcvb(nproc*arrsize) )
end if

! zero arrays
flag=0
mbnum=0
nam=0

! loop over all beads on each node 
do i=1,natms
    i1=corr(i,1)
    nb=cn(i1,0)

    ! loop over all bonds of the bead
    do k=1,nb
        i2=cn(i1,k)
        arrtype=bcorrt(i2)

        ! check if there is information about bonded bead
        if (arrtype/=0) then
            j=bcorr(i2)

            ! calculate distance between beads
            rxij = rx(i,1) - rx(j,arrtype)
            ryij = ry(i,1) - ry(j,arrtype)
            rzij = rz(i,1) - rz(j,arrtype)

            ! apply PBC
            if (rxij>dlxa/2) then 
                rxij=rxij-dlxa
            elseif (rxij<-dlxa/2) then
                rxij=rxij+dlxa
            end if
            if (ryij>dlya/2) then 
                ryij=ryij-dlya
            elseif (ryij<-dlya/2) then
                ryij=ryij+dlya
            end if
            if (rzij>dlza/2) then 
                rzij=rzij-dlza
            elseif (rzij<-dlza/2) then
                rzij=rzij+dlza
            end if
            rijsq=sqrt(rxij**2+ryij**2+rzij**2)

            ! calculate forces
            fx(i) = fx(i) - k_bond(kindpall(i1),kindpall(i2)) * rxij*(1-k_eq(kindpall(i1),kindpall(i2))/rijsq)
            fy(i) = fy(i) - k_bond(kindpall(i1),kindpall(i2)) * ryij*(1-k_eq(kindpall(i1),kindpall(i2))/rijsq)
            fz(i) = fz(i) - k_bond(kindpall(i1),kindpall(i2)) * rzij*(1-k_eq(kindpall(i1),kindpall(i2))/rijsq)

            ! calculate virial for pressure calculation
            press(1)=press(1)-k_bond(kindpall(i1),kindpall(i2)) * rxij*(1-k_eq(kindpall(i1),kindpall(i2))/rijsq)*rxij*0.5
            press(2)=press(2)-k_bond(kindpall(i1),kindpall(i2)) * ryij*(1-k_eq(kindpall(i1),kindpall(i2))/rijsq)*ryij*0.5
            press(3)=press(3)-k_bond(kindpall(i1),kindpall(i2)) * rzij*(1-k_eq(kindpall(i1),kindpall(i2))/rijsq)*rzij*0.5

        ! if there is no information about bonded bead - put it to the list of lost connections
        else
            flag=1
            mbnum=mbnum+1
            mb(mbnum,1)=i2
            mb(mbnum,2)=i
            nam=nam+1
            if (nam>arrsize) call error(160)
            sndb(nam)=real(i2)
        end if
    end do
end do

! gather number of lost connections
call MPI_Allgather (nam, 1, MPI_integer, rcount, 1, MPI_integer, MPI_Comm_world,ierr )
displ(1)=0
do i=2,nproc
    displ(i)=displ(i-1)+rcount(i-1)
end do
allnum=displ(nproc)+rcount(nproc)

! gather list of lost connections
call MPI_Allgatherv(sndb, nam,MPI_real4, rcvb, rcount, displ, MPI_real4, MPI_Comm_world,ierr)
nam=0

! look for information about lost beads on each node
do k=1,allnum
    if (bcorrt(int(rcvb(k)))==1) then
        i=bcorr(int(rcvb(k)))
        nam=nam+4
        if (nam>arrsize) call error(160)
        sndb(nam-3)=rcvb(k)
        sndb(nam-2)=rx(i,1)+(ncol-1)*dlx
        sndb(nam-1)=ry(i,1)+(nrow-1)*dly
        sndb(nam)=rz(i,1)
    end if
end do    

! gather number of found connections 
call MPI_Allgather (nam, 1, MPI_integer, rcount, 1, MPI_integer, MPI_Comm_world,ierr )
displ(1)=0
do i=2,nproc
    displ(i)=displ(i-1)+rcount(i-1)
end do
allnum=displ(nproc)+rcount(nproc)

! gather linformation about found connections
call MPI_Allgatherv(sndb, nam,MPI_real4, rcvb, rcount, displ, MPI_real4, MPI_Comm_world,ierr)

! look for information about lost beads in found connections
! loop over all found connections
do i=0,allnum-1,4
    i2b=rcvb(i+1)
    shft=0

    ! loop over all lost beads
    do j=1,mbnum
        if (mb(j-shft,1)==i2b) then
            i1=mb(j-shft,2)
            i3=corr(i1,1)

            ! calculate distance between beads
            rxij = rx(i1,1) - rcvb(i+2)+(ncol-1)*dlx
            ryij = ry(i1,1) - rcvb(i+3)+(nrow-1)*dly
            rzij = rz(i1,1) - rcvb(i+4)

            ! apply PBC
            if (rxij>dlxa/2) then 
                rxij=rxij-dlxa
            elseif (rxij<-dlxa/2) then
                rxij=rxij+dlxa
            end if
            if (ryij>dlya/2) then 
                ryij=ryij-dlya
            elseif (ryij<-dlya/2) then
                ryij=ryij+dlya
            end if
            if (rzij>dlza/2) then 
                rzij=rzij-dlza
            elseif (rzij<-dlza/2) then
                rzij=rzij+dlza
            end if
            rijsq=sqrt(rxij**2+ryij**2+rzij**2)

            ! calculate forces
            fx(i1) = fx(i1) - k_bond(kindpall(i3),kindpall(i2b)) * rxij*(1-k_eq(kindpall(i3),kindpall(i2b))/rijsq)
            fy(i1) = fy(i1) - k_bond(kindpall(i3),kindpall(i2b)) * ryij*(1-k_eq(kindpall(i3),kindpall(i2b))/rijsq)
            fz(i1) = fz(i1) - k_bond(kindpall(i3),kindpall(i2b)) * rzij*(1-k_eq(kindpall(i3),kindpall(i2b))/rijsq)

            ! calculate virial for pressure calculation
            press(1)=press(1)-k_bond(kindpall(i3),kindpall(i2b)) * rxij*(1-k_eq(kindpall(i3),kindpall(i2b))/rijsq)*rxij*0.5
            press(2)=press(2)-k_bond(kindpall(i3),kindpall(i2b)) * ryij*(1-k_eq(kindpall(i3),kindpall(i2b))/rijsq)*ryij*0.5
            press(3)=press(3)-k_bond(kindpall(i3),kindpall(i2b)) * rzij*(1-k_eq(kindpall(i3),kindpall(i2b))/rijsq)*rzij*0.5
            mb(j-shft,:)=mb(mbnum-shft,:)
            shft=shft+1
        end if
    end do
    mbnum=mbnum-shft
end do
if (mbnum/=0) then
    call error(161)
end if

end

subroutine vcontrol()
! #################################################################
! #                                                               #
! #    subroutine 17:                                             #
! #    prevent moving of the center of mass                       #
! #                                                               #
! #################################################################
implicit none
common /vxyz/ vx,vy,vz
common /na/ natms,natmsfan
common /npall/ npall
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
real*4 vom(3)
integer*4 i,natms,natmsfan,npall

! velocity of the center of mass
vom=0
do i =1, natms
    vom(1) = vom(1) + vx(i)
    vom(2) = vom(2) + vy(i)
    vom(3) = vom(3) + vz(i)
end do
call gsum(vom,3)
vom = vom/real(npall)

! stop the center of mass
do i =1, natms
    vx(i)=vx(i)-vom(1)
    vy(i)=vy(i)-vom(2)
    vz(i)=vz(i)-vom(3)
end do

end

subroutine writerst()
! #################################################################
! #                                                               #
! #    subroutine 18:                                             #
! #    write a restart file                                       #
! #                                                               #
! #################################################################
implicit none
save

common /cell/ dlx,dly,dlz,rcut,kpbc
common /cella/ dlxa,dlya,dlza
common /npall/ npall
common /na/ natms,natmsfan
common /mpidpd/ rank, nproc
common /rxyz/ rx,ry,rz
common /vxyz/ vx,vy,vz
common /ptype/ kindp
common /rst/ nbond
common /bnd/ cn 
common /correl/ corr,bcorr,bcorrt
common /top/ cln,nrow,ncol 
common /vale/ val 
common /den/ rho  
common /ptypeall/ kindpall

real*4 dlx,dly,dlz,rcut
real*4 dlxa,dlya,dlza
integer*4 npall,kpbc,natms,natmsfan
integer*4 rank,nproc,cln,nrow,ncol
integer*4 i,j,ip,nre,ipos
integer*4 nbond,nangles
real*4 rho
character*16 cre
integer b
logical init/.true./

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:)
real*4,pointer :: vx(:),vy(:),vz(:) 
integer*1,pointer :: kindp(:,:)
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:) 
integer*4,pointer :: cn(:,:)
integer*1,pointer :: val(:)
integer*1,pointer :: kindpall(:)

! set initial values
if (init) then
    init=.false.
    nre=0
    nangles=0
end if

! calculate current restart file number
nre=nre+1
write(cre,'(i16)') nre
do ipos = 1,16
    if (cre(ipos:ipos)/=' ') exit
end do

! write system size and number of particles
if (rank==0) then
    open( 11, file = 'restart'//cre(ipos:16)//'.dat' )
    write(11,'(i8,f8.4)') npall,rho
    write(11,'(3f20.12)') dlxa,dlya,dlza
    close(11, status = 'keep') 
end if

! write information about beads: coordinates, types, valencies
do ip=0,nproc-1
    if (rank==ip) then
        open( 11, file = 'restart'//cre(ipos:16)//'.dat',position= 'append' )
        do i = 1, natms
            write(11,'(i8,2i4,3f14.6)')corr(i,1),val(corr(i,1)),kindp(i,1),rx(i,1)+(ncol-1)*dlx,ry(i,1)+(nrow-1)*dly,rz(i,1)
        end do
        close(11,status = 'keep')
    end if
    call barrier
end do

if (rank==0) then
    open( 11, file = 'restart'//cre(ipos:16)//'.dat',position= 'append' )

    ! write bonds
    write(11,'(a9,i9)') ' bonds:  ',nbond
    do i = 1, npall
        do j=1,cn(i,0)
            if(i<cn(i,j)) write(11,'(2i8)') i, cn(i,j)
        end do    
    end do

        ! write angles
    write(11,'(a9,i9)') ' angles: ',nangles
    close(11, status = 'keep')
end if

end


subroutine pressure()
! #################################################################
! #                                                               #
! #    subroutine 24:                                             #
! #    calculate pressure using virial theorem                    #
! #                                                               #
! #################################################################
implicit none
include "mpif.h"

common /vxyz/ vx,vy,vz
common /cella/ dlxa,dlya,dlza
common /na/ natms,natmsfan
common /prs/ press

real*4,pointer :: vx(:),vy(:),vz(:) 

real*4 volm
integer*4 i
real*4 press(3)
real*4 dlxa,dlya,dlza
integer*4 natms,natmsfan
integer ierr

! cell volume
volm=dlxa*dlya*dlza

! calculate kinetic contribution to the virial
do i=1,natms
    press(1)=press(1)+vx(i)*vx(i)
    press(2)=press(2)+vy(i)*vy(i)
    press(3)=press(3)+vz(i)*vz(i)
end do

! collect contributions of all nodes
call mpi_allreduce(mpi_in_place,press,3,mpi_real,mpi_sum,mpi_comm_world,ierr)

! calculate pressure
press=press/volm

end

Subroutine distij(i,j,arrtype,rij)
! #################################################################
! #                                                               #
! #    subroutine 28:                                             #
! #    calculate distance between particles                       #
! #                                                               #
! #################################################################
implicit none

common /cella/ dlxa,dlya,dlza
common /rxyz/ rx,ry,rz 

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:)

real*4 dlxa,dlya,dlza,rxb,ryb,rzb,rij
integer*4 arrtype,i,j

rxb=rx(i,1)-rx(j,arrtype)
ryb=ry(i,1)-ry(j,arrtype)
rzb=rz(i,1)-rz(j,arrtype)

! apply PBC
rxb=rxb-dlxa*nint(rxb/dlxa)
ryb=ryb-dlya*nint(ryb/dlya)
rzb=rzb-dlza*nint(rzb/dlza)

rij=rxb**2+ryb**2+rzb**2

end


subroutine sort(nneib,reaclist,reaclistc,reaclistart,distlist)
! #################################################################
! #                                                               #
! #    subroutine 29:                                             #
! #    sort list of possible connections(bubble sort)             #
! #                                                               #
! #################################################################   
implicit none
      
integer*4 nneib,kl,flag,i
integer*4 reaclist(1000),reaclistc(1000),reaclistart(1000)
real*4 distlist(1000),lk

! classical bubble sort      
334   flag=0
do i=1,nneib-1
    if (distlist(i)>distlist(i+1)) then
        flag=1
        lk=distlist(i)
        distlist(i)=distlist(i+1)
        distlist(i+1)=lk
        
        kl=reaclist(i)
        reaclist(i)=reaclist(i+1)
        reaclist(i+1)=kl
        
        kl=reaclistc(i)
        reaclistc(i)=reaclistc(i+1)
        reaclistc(i+1)=kl
        
        kl=reaclistart(i)
        reaclistart(i)=reaclistart(i+1)
        reaclistart(i+1)=kl
    end if
end do
if (flag==1) go to 334

end

subroutine rconf(steps1,steps2,qstp)
! #################################################################
! #                                                               #
! #    subroutine 32:                                             #
! #    read input script                                          #
! #                                                               #
! #################################################################  
implicit none

common /cella/ dlxa,dlya,dlza 
common /den/ rho 
common /datatech/ alpha,sigma,gamma
common /time/ dt
common /shwi/ shwi
common /mpidpd/ rank, nproc 
common /bondf/ k_eq,k_bond 
common /bmult/ mult
common /lambd/ lambda

real*4,pointer :: alpha(:,:)
character*2,pointer :: bead_type(:)
real*4,pointer :: k_eq(:,:),k_bond(:,:)

character*256 str
character*2 el

integer*4 i,j,k
integer*1 rtype,snum,spnum,s1,s2,atempl
real*4 dlxa,dlya,dlza,dt,sigma,gamma
integer*4 steps1,steps2
real*4 alph,prob,probr
real*4 rho,shwi
integer*4 qstp
integer*4 rank,nproc,mult
logical exists
real*4 lambda

! set default values
! * technical values
sigma=3.0
shwi=1.0
steps1=1000
steps2=100000
spnum=0
mult=5
dt=0.04
qstp=1000
atempl=0

! allocate all arrays
allocate(k_eq(10,10),k_bond(10,10),alpha(10,10),bead_type(10))
alpha=25.0
k_eq=0.5
k_bond=40.0
lambda=0.65

alpha(1,2)=26.75
alpha(2,1)=26.75
rho=3.0
bead_type(1)='C'
bead_type(2)='O'


gamma=0.50*sigma*sigma

! construct initial system
!if (rtype==1.and.rank==0) call construct(fnames,conc,place,snum,anglestempl,atempl,vallist,bead_type,ntlength,sphrad,graftdens)
if (rank==0) call read_geo()
call barrier


return
1 call error(323)    !stop 'error:spnum=0'
4 call error(326)    !    stop 'error:n_species error'
5 call error(327)    !    stop 'error:box size error'
6 call error(328)    !    stop 'error:prname error'
7 call error(329)    !    stop 'error:dt error'
8 call error(3210)    !    stop 'error:t1 error'
9 call error(3211)    !    stop 'error:t1 error'
10 call error(3212)    !   stop 'error:t1 error'
11 call error(3213)    !   stop 'error:tch error'
12 call error(3214)    !   stop 'error:alpha error'
14 call error(3216)    !   stop 'there is no dpdconf.dat file! '
15 call error(3217)    !   stop 'error:density error!'
16 call error(3218)    !   stop 'output_freq error'
18 call error(3220)    !   stop 'btempl error'
19 call error(3221)    !   stop 'val error'
23 call error(3225)    !   stop 'restart error'
25 call error(3227)    !   stop 'bead_type error'
26 call error(3228)    !   stop 'shwi error'
27 call error(3229)    !   stop 'shwi error'
28 call error(3230)    !   stop 'shwi error'
29 call error(3231)
30 call error(3232)
end

function trimtext(string)
! #################################################################
! #                                                               #
! #    function 33:                                               #
! #    find last non-blank character                              #
! #                                                               #
! #################################################################  
implicit none

integer i,size,len
integer last,trimtext
character*1 null
character*(*) string

! move forward through the string, one character
! at a time, looking for first null character

trimtext = 0
null = char(0)
size = len(string)
last = size
do i = 1, size
    if (string(i:i) == null) then
        last = i - 1
        goto 10
    end if
end do
10 continue

! move backward through the string, one character
! at a time, looking for first non-blank character

do i = last, 1, -1
    if (string(i:i) > ' ') then
        trimtext = i
        goto 20
    end if
end do

20 continue

end


subroutine rndplace( nat,rxt,ryt,rzt,cn,dlxa,dlya,dlza, place )
! #################################################################
! #                                                               #
! #    subroutine 36:                                             #
! #    place a single DPD molecule in the periodic box            #
! #                                                               #
! ################################################################# 
implicit none      

real*4 rxt(*),ryt(*),rzt(*)
integer*4 cn(nat,0:*)
integer*4 nat, place
real*4 dlxa,dlya,dlza,uni
integer*1 id(nat)
integer*4 i,k,j
real*4 dx,dy,dz,r,scale

! generate initial beads coordinates 
do i = 1, nat
    id(i)=1
    rxt(i)=dlxa*(uni())
    ryt(i)=dlya*(uni())
    rzt(i)=dlza*(uni())
end do

! correct initial coordinates considering connection matrix
do i = 1, nat-1

    ! exit if the bead has no bonds
    if ( cn(i,0) == 0 ) then
        id(i)=-1
        cycle
    end if

    ! place all connected to the ith bead beads next to it
    !do k = 1, 2!cn(i,0)

        j = i+1
        ! don't move jth bead if it has already been moved 
        if (id(j)==-1) cycle
        dx = ( uni() - 0.5 )
        dy = ( uni() - 0.5 )
        dz = ( uni() - 0.5 )
        r = sqrt( dx**2 + dy**2 + dz**2 )
        if ( r==0.0 ) then
		r = 1.0
		end if
		scale = 1.0 / r
        dx = dx * scale
        dy = dy * scale
        dz = dz * scale

        if (rxt(i)+dx < dlxa .and. rxt(i)+dx > 0) then
            rxt(j)=rxt(i)+dx
        else
            rxt(j)=rxt(i)-dx
        end if

        if (ryt(i) + dy < dlya .and. ryt(i) + dy > 0) then
            ryt(j) = ryt(i) + dy
        else
            ryt(j) = ryt(i) - dy
        end if

        if (rzt(i) + dz < dlza .and. rzt(i) + dz > 0) then
            rzt(j) = rzt(i) + dz
        else
            rzt(j) = rzt(i) - dz
        end if

		!write(*,*) dx, dy, dz, scale, r, place
        id(j) = -1
    !end do
end do

end

subroutine error(flag)
! #################################################################
! #                                                               #
! #    subroutine 37:                                             #
! #    stop the program due to error                              #
! #                                                               #
! ################################################################# 
implicit none
common /mpidpd/ rank, nproc

integer*4 flag
integer*4 rank,nproc

if( flag==60) then
    stop 'error 60:a bead has too high speed along X'
end if
if( flag==61) then
    stop  'error 61:a bead has too high speed along Y'
end if
if( flag==150) then
    call barrier
    if(rank==0)stop 'error 150:domain size along X is less than Rcut'
    stop
end if
if( flag==151) then
    call barrier
    if(rank==0)stop 'error 151:domain size along Y is less than Rcut'
    stop
end if
if( flag==152) then
    call barrier
    if(rank==0)stop 'error 152:a bead has been lost during reading restart file'
    stop
end if
if( flag==153) then
    call barrier
    if(rank==0)stop 'error 153:restart.dat file does not exist'
    stop
end if
if( flag==154) then
    stop 'error 154:error during reading restart.dat file'
end if
if( flag==160) then
    stop 'error 160:communication array overflow in forcesbond subroutine'
end if
if( flag==320) then
    stop 'error 320: input script error, unrecognized run_type'
end if
if( flag==321) then
    stop 'error 321: input script error, run_type should be specified before all other commands'
end if
if( flag==322) then
    stop 'error 322: input script error, number of species should be specified before input_file command'
end if
if( flag==327) then
    stop 'error 327: input script error, box_size command error'
end if
if( flag==329) then
    stop 'error 329: input script error, timestep command error'
end if
if( flag==3210) then
    stop 'error 3210: input script error, steps_1 command error'
end if
if( flag==3211) then
    stop 'error 3211: input script error, steps_2 command error'
end if
if( flag==3214) then
    stop 'error 3214: input script error, pair_coeff command error'
end if
if( flag==3217) then
    stop 'error 3217: input script error, density command error'
end if
if( flag==3218) then
    stop 'error 3218: input script error, output_freq command error'
end if
if( flag==3220) then
    stop 'error 3220: input script error, bond_template command error'
end if
if( flag==3225) then
    stop 'error 3225: input script error, restart command error'
end if
if( flag==3227) then
    stop 'error 3227: input script error, bead_type command error'
end if
if( flag==3228) then
    stop 'error 3228: input script error, shwi command error'
end if
if( flag==340) then
    stop 'error 340:element type is not set'
end if
write(*,*) 'error with flag', flag
stop 'unspecified error'
end


subroutine gsum(arr,n)
! #################################################################
! #                                                               #
! #    subroutine 38:                                             #
! #    shortcut to mpi_allreduce                                  #
! #                                                               #
! ################################################################# 
implicit none
include 'mpif.h'
     
real*4 arr(n),arr1(n)
integer*4 n,i
integer ierr
     
call MPI_allreduce(arr,arr1,n,MPI_real,MPI_SUM,MPI_COMM_WORLD,ierr)
do i=1,n
    arr(i)=arr1(i)
end do

end

subroutine barrier()
! #################################################################
! #                                                               #
! #    subroutine 39:                                             #
! #    shortcut to mpi_barrier                                    #
! #                                                               #
! #################################################################
include "mpif.h"
integer ierr

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end
