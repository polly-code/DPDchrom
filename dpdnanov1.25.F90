
implicit none
include "mpif.h"

common /cell/ dlx,dly,dlz,rcut,kpbc !box size on each node (dlx,dly,dlz), cut radius (rcut), type of PBC(kpbc)
common /cella/ dlxa,dlya,dlza !complete cox size
common /sizes/ sz !size of main arrays (like r and v)
common /bxsize/ ilx,ily,ilz !number of box bins 
common /na/ natms,natmsfan !current number of "real" particles(natms), current number of right "ghost" particles
common /nab/ natmsfanbond !current number of left "ghost" particles. necessary for bond forces calculation
common /data/ alpha,sigma,gamma !dpd parameters 
common /time/ dt ! timestep
common /comm/ ssz,fsz !size of communication arrays, ssz is for moving atoms between nodes, fsz - for "ghost" atoms 
common /npall/ npall !total number of atoms in the system
common /shwi/ shwi !transfer radius for "ghost" atoms 
common /correl/ corr,bcorr,bcorrt !arrays that store global atom numbers. corr converts atom number on the node to the global number, bcorr does the opposite, bcorrt shows atom state ("real,"ghost" or on the other node)
common /den/ rho !density 
common /bondf/ k_eq,k_bond !bond parameters - eq.length, stiffness
common /anglef/ k_eqa,k_angle !angle parameters - eq.value, energy coefficient
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
common /rst/ nbond,nangles !number of bonds,number of angles
common /typel/ typelist,ntype !type of the atoms that present in the system and number of such types
common /ptypeall/ kindpall !particle types
common /top/ cln,nrow,ncol ! number of nodes in 1 row,row and column number of node
common /angl/ ca,angles ! main angle arrays
common /prs/ press ! pressure
common /vale/ val !valencies
common /probcb/ probc,probb !probabilities of chemical reactions
common /ch/ rchem !radius of chemical reactions
common /bmult/ mult !angles/bond snd rcv arrays size multiplier
common /defor/ nave ! number of interations to average stress-strain curve
common /lambd/ lambda

real*4,pointer :: sndr(:,:),sndl(:,:),fanl(:,:),fanr(:,:),rcv(:),rcvf(:) !size ssz,4;ssz,4;fsz,4;fsz,4;ssz*8,fsz*8
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) ! size sz,3
real*4,pointer :: vwx(:,:),vwy(:,:),vwz(:,:) ! size sz,2
real*4,pointer :: fx(:),fy(:),fz(:) ! size sz
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
integer*1,pointer :: kindp(:,:) ! size sz,3
integer*4,pointer :: corr(:,:),bcorr(:)! size sz,3 and npall
integer*1,pointer :: bcorrt(:) !size npall
integer*4,pointer :: fann(:,:) !size fsz/8,4
real*4,pointer :: fsend(:) ! size 3*fsz/2
integer*4,pointer :: cn(:,:) !size npall,0:maximum number of bonds
integer*1,pointer :: kindpall(:) ! size npall
integer*4,pointer :: ca(:,:) !size npall, 0: maximum number of angles
integer*4,pointer :: angles(:) !size 3*nangles
integer*1,pointer :: val(:) !size npall
real*4,pointer :: alpha(:,:),probc(:,:),probb(:,:) !size 10,10
real*4,pointer :: k_eq(:,:),k_bond(:,:) !size 10,10
real*4,pointer :: k_eqa(:,:,:),k_angle(:,:,:) !size 10,10,10

integer*4 nbond,nangles
integer*4 rank,nproc,cln,nrow,ncol
real*4 uni
integer*4 sz,ssz,fsz,npall
integer*4 ilx,ily,ilz
integer ierr
integer*4 steps1,steps2,steps3,qstp,rststp,sqfstp,vtkstp,nstp,chemstp,stepsdef,ndef,nave
integer*4 kpbc,rst,sqf,vtk
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
real*4 press(3),tau,tai
real*4 rchem
real*4 lambda
integer*1 def_res

! MPI initialization
call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)

! randon number generator initialization
vl=uni()

! read input script
call rconf(steps1,steps2,steps3,qstp,rst,rststp,vtk,vtkstp,sqf,sqfstp,chemstp,stepsdef,ndef,tau,tai,def_res)

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
    write(*,*) 'dpdnanov1.25k'
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
    write(*,*) 'Number of steps in 3rd stage :',steps3
    write(*,*) 'Number of steps to deform    :',stepsdef
    write(*,*)
    write(*,*) 'DPD Parameters'
    write(*,*) 'Sigma                        :',sigma
    write(*,*) 'Gamma                        :',gamma
    do i=1,ntype
        do j=i,ntype
            write(*,*) 'DPD coeff btw',typelist(i),'and',typelist(j),'is',alpha(typelist(i),typelist(j))
        end do
    end do

    do i=1,ntype
        do j=i,ntype
            write(*,*) 'Chem probabilities btw',typelist(i),'and',typelist(j),'are',probc(typelist(i),typelist(j)),probb(typelist(i),typelist(j))
        end do
    end do

    write(*,*)
end if

! output vtk file
if(vtk==1) call dpd2vtk()

! calculate forces
fsend=0
call forces()
call forcesbond()
call forcesangle()
if (rank==0) start=mpi_wtime()

! start 1st stage with soft potentials
alpha=alpha/5.0
k_bond=k_bond/5.0
do nstp=1,steps1
    
    ! integrate beads positions and velocities
    call intr()

    ! calculate forces
    call forces()
    call forcesbond()
    call forcesangle()

    ! integrate beads velocoties
    call intv()

    ! output vtk file
    if (mod(nstp,vtkstp)==0.and.vtk==1) call dpd2vtk()
    
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
        if (rank==0) write(*,*)'stage 1', nstp,vl,maxval(abs(vall))
    end if
    
    ! control  CM velocily
    if (mod(nstp,100)==0) call vcontrol()

    ! output restart file
    if (mod(nstp,rststp)==0.and.rst==1) call writerst()
end do

! finish 1st pahse
alpha=alpha*5.0
k_bond=k_bond*5.0

! output time
if (rank==0) write(*,*) 'stage 1 took',mpi_wtime()-start,'s'

! output VTK file
if (vtk==1) call dpd2vtk()

! output restart file
if (rst==1) call writerst()

! start 2nd stage - chem stage
if (rank==0) start=mpi_wtime()
do nstp=1,steps2

    ! integrate beads positions and velocities
    call intr()

    ! calculate forces
    call forces()
    call forcesbond()
    call forcesangle()

    ! integrate beads velocoties
    call intv()

    ! output VTK file
    if (mod(nstp,vtkstp)==0.and.vtk==1) call dpd2vtk()

    ! perform reactions
    if (mod(nstp,chemstp)==0) call chem()

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

    ! control  CM velocily
    if (mod(nstp,100)==0) call vcontrol()

    ! output restart file
    if (mod(nstp,rststp)==0.and.rst==1) call writerst()
end do

! finish 2nd stage 
! output time
if (rank==0) write(*,*) 'stage 2 took',mpi_wtime()-start,'s'

! output VTK file
if (vtk==1) call dpd2vtk()

! output restart file
if (rst==1) call writerst()

! start 3rd stage - equilibration and sq stage
if (rank==0) start=mpi_wtime()
do nstp=1,steps3

    ! integrate beads positions and velocities
    call intr()

    ! calculate forces
    call forces()
    call forcesbond()
    call forcesangle()

    ! integrate beads velocoties
    call intv()

    ! output VTK file
    if (mod(nstp,vtkstp)==0.and.vtk==1) call dpd2vtk()

    ! calculate structure factor
    if (mod(nstp,sqfstp)==0) then
        if (sqf==1) call sqs()
    end if

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
        if (rank==0) write(*,*)'stage 3', nstp,vl,maxval(abs(vall)),sum(press)/3.0
    end if

    ! control  CM velocily
    if (mod(nstp,100)==0) call vcontrol()

    ! output restart file
    if (mod(nstp,rststp)==0.and.rst==1) call writerst()
end do

! finish 3rd stage 
! output time
if (rank==0) write(*,*) 'stage 3 took',mpi_wtime()-start,'s'

! output VTK file
if (vtk==1) call dpd2vtk()

! output restart file
if (rst==1) call writerst()

! start deformation
if (rank==0) start=mpi_wtime()
dlxtar=dlxa+def_res*tau

! loop over number of deformed states
do curdef=0,ndef
    do nstp=1,stepsdef
        
        ! integrate beads positions and velocities
        call intr()

        ! calculate forces
        call forces()
        call forcesbond()
        call forcesangle()

        ! integrate velocities
        call intv()

        ! output VTK file
        if (mod(nstp,vtkstp)==0.and.vtk==1) call dpd2vtk()

        ! perform system affine deformation
        if(tau>0.and.dlxa<dlxtar.and.mod(nstp,100)==0) call deform(tau,tai)
        if(tau<0.and.dlxa>dlxtar.and.mod(nstp,100)==0) call deform(tau,tai)
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
            if(rank==0.and.nstp>(stepsdef-nave*qstp))call pressout(0)
            if (rank==0) write(*,*)'stage def', nstp,vl,maxval(abs(vall)),press(1)
        end if

        ! control  CM velocily
        if (mod(nstp,100)==0) call vcontrol()

        ! output restart file
        if (mod(nstp,rststp)==0.and.rst==1) call writerst()

    end do

    ! output Stress-Strain curve
    if(rank==0)call pressout(1)

    ! output restart file
    if (rst==1) call writerst()

    ! move to the next deformation degree
    dlxtar=dlxtar+tau
end do

! output time
if (rank==0) write(*,*) 'stage def took',mpi_wtime()-start,'s'

! output VTK file
if (vtk==1) call dpd2vtk()

! delete arrays
call barrier
nullify(rx,ry,rz,vwx,vwy,vwz,fx,fy,fz,vx,vy,vz,kindp,corr,bcorr,bcorrt,sndr,sndl,fanl,rcv,rcvf,fann,fsend)

! end mpi
call mpi_finalize(ierr)

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

common /data/ alpha,sigma,gamma
common /time/ dt
common /sizes/ sz
common /fxyz/ fx,fy,fz
common /vwxyz/ vwx,vwy,vwz
common /ptype/ kindp
common /commf/ fsend,rcou
common /prs/ press 

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

! set initial values
if (init) then
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

! set initial values
if (init) then
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

subroutine dpd2vtk()
! #####################################################################
! #                                                                   #
! #   subroutine 12:                                                  #
! #   output VTK files                                                #
! #                                                                   #
! #####################################################################
implicit none
save
include "mpif.h"

common /ptype/ kindp
common /na/ natms,natmsfan
common /cell/ dlx,dly,dlz,rcut,kpbc
common /cella/ dlxa,dlya,dlza
common /bxsize/ ilx,ily,ilz
common /rxyz/ rx,ry,rz
common /mpidpd/ rank, nproc
common /den/ rho
common /typel/ typelist,ntype
common /top/ cln,nrow,ncol

logical init/.true./
integer*4 mw,i,j,k,ijk,m
real*4 sw,rcut
integer*4   iw(0:26), jw(0:26), kw(0:26)
real*4      w(0:26)
integer*4 b,kpbc
integer*1 typelist(10),ntype,itype
integer*4 natms,natmsfan
integer*4 frame
character*4 fname(10000)
character*1 iat(10)
character*1 number(10)
integer*4 ni,nj
real*4 dlx,dly,dlz
real*4 dlxa,dlya,dlza
integer*4 ilx,ily,ilz
real*4 x0,y0,z0
integer*4 ix,iy,iz
integer*4 rank, nproc,cln,nrow,ncol
real*4 di,d(int(dlxa))
real*4 rho
integer ierr

integer*1,pointer :: kindp(:,:) 
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 

integer*4   nt(int(dlxa),int(dlya),ilz), nta(int(dlxa),int(dlya),ilz)
integer*4   ic(0:int(dlxa)+1), jc(0:int(dlya)+1), kc(0:ilz+1)
data number / '0','1','2','3','4','5','6','7','8','9' /
data iat/'1','2','3','4','5','6','7','8','9','0'/

! set intial values   
if ( init ) then
    init = .false.
    mw = 0
    sw = 1.0
    iw(0) = 0
    jw(0) = 0
    kw(0) = 0
    w(0)  = 1.0

    ! set weighting coefficients
    do i = -1, 1
        do j = -1, 1
            do k = -1, 1
                ijk = i*i + j*j + k*k
                if ( ijk > 0 ) then
                    mw = mw + 1
                    iw(mw) = i
                    jw(mw) = j
                    kw(mw) = k
                    if ( ijk == 1 ) w(mw) = 1.d0/( 1.+sqrt( 1. ) )
                    if ( ijk == 2 ) w(mw) = 1.d0/( 1.+sqrt( 2. ) )
                    if ( ijk == 3 ) w(mw) = 1.d0/( 1.+sqrt( 3. ) )
                    sw = sw + w(mw)
                end if
            end do
        end do
    end do

    !convert frame number to character
    frame = 0
    do i = 1, 10
        do j = 1, 10
            do k = 1, 10
                do m=1, 10
                    frame = frame + 1
                    fname(frame) = number(i)//number(j)//number(k)//number(m)
                end do
            end do
        end do
    end do
    frame = 0
end if
frame = frame + 1

! exit program if there are too many frames
if ( frame > 9999 ) call error(120)

! create subcell map considering boundary conditions
ni =  int( dlxa )
nj =  int( dlya )
do i = 0, ni + 1
    ic(i) = mod( i, ni )
    if ( ic(i) == 0 ) ic(i) = ni
end do
do j = 0, nj + 1
    jc(j) = mod( j, nj )
    if ( jc(j) == 0 ) jc(j) = nj
end do
do k = 0, ilz + 1
    kc(k) = mod( k, ilz )
    if ( kc(k) == 0 ) kc(k) = ilz
end do

! loop over all type of particles in the simulation cell
do itype = 1, ntype
    ijk = typelist(itype)
    nt=0
    nta=0

    ! each node calculates the number of beads of chosen type in subcells
    do i = 1, natms
        if ( kindp(i,1) /= ijk ) cycle
        x0 = rx(i,1)+dlx*(ncol-1)
        y0 = ry(i,1)+dly*(nrow-1)
        z0 = rz(i,1)
        ix=min(int(x0),int(dlxa)-1)+1
        iy=min(int(y0),int(dlya)-1)+1
        iz=min(int(z0),ilz-1)+1
        nt(ix,iy,iz) = nt(ix,iy,iz) + 1
    end do

    ! collect number of beads in all subcells on the node 0 
    call MPI_REDUCE(nt, nta, ni*nj*ilz, mpi_integer, mpi_sum, 0, mpi_comm_world,ierr)

    ! output vtk file considering weighting coefficients
    if (rank==0) then
        open(13, file = iat(typelist(itype))//'-'//fname(frame)//'.vtk',status='replace')
        write( 13, '(a)' )      '# vtk DataFile Version 3.0'
        write( 13, '(a)' )      'Uni-Ulm Simuation Group: rgba'
        write( 13, '(a)' )      'ASCII'
        write( 13, '(a)' )      'DATASET STRUCTURED_POINTS'
        write( 13, '(a,3i4)' )  'DIMENSIONS', ni,nj,ilz
        write( 13, '(a)' )      'ASPECT_RATIO 1 1 1'
        write( 13, '(a)' )      'ORIGIN 0 0 0'
        write( 13, '(a,i8)' )   'POINT_DATA', ni*nj*ilz
        write( 13, '(a)' )      'SCALARS volume_scalars int 1'
        write( 13, '(a)' )      'LOOKUP_TABLE default'
        do k = 1, ilz
            do j = 1, nj
                do i = 1, ni
                    di = 0.0

                    ! calculate final subcell beads density considering contribution of neighbouring subcells
                    do m = 0, 26
                        di = di +nta(ic(i+iw(m)),jc(j+jw(m)),kc(k+kw(m)))*w(m)
                    end do
                    d(i) = di / sw / rho
                end do
                write(13,'(128i4)') ( int( 500*d(i) ), i = 1, ni )
            end do
        end do
        close(13, status = 'keep')
    end if
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
real*4,pointer :: fx(:),fy(:),fz(:) ! size sz
real*4,pointer :: vx(:),vy(:),vz(:) ! size sz 
real*4 dt,dt2
integer*4 i,natms,natmsfan
logical init/.true./

! set initial value
if (init) then
    dt2 = 0.50*dt
    init=.false.
end if
      
! integrate velocities
do i=1, natms
    vx(i) = vx(i) + dt2*fx(i)
    vy(i) = vy(i) + dt2*fy(i)
    vz(i) = vz(i) + dt2*fz(i)
end do

end

subroutine sqs()
! #################################################################
! #                                                               #
! #    subroutine 14:                                             #
! #    calculate static structure factor                          #
! #                                                               #
! #################################################################
implicit none
include 'mpif.h'
save

common /cella/ dlxa,dlya,dlza
common /cell/ dlx,dly,dlz,rcut,kpbc
common /na/ natms,natmsfan
common /rxyz/ rx,ry,rz
common /ptype/ kindp
common /mpidpd/ rank, nproc
common /top/ cln,nrow,ncol

integer kmax,j,numvect,i,k,len,i2st,j2st,k2st,atom,counter,i2,j2,k2,i4,l
parameter (kmax=33)
parameter (numvect=2000)
logical init /.true./
real*4 twopi,sq(numvect,10)
real*4 minvec1,minvec2,minvec3,rn1,rn2,rn3
real*4 dlxa,dlya,dlza
real*4 dlx,dly,dlz,rcut
integer*4 kpbc
real*4 rxb,ryb,rzb
integer*4 natms,natmsfan
integer*4 rank,nproc,cln,nrow,ncol
integer num(numvect),nat(10),nata(10),ierr,kinda
complex expsum(0:kmax,0:kmax,8,10),expsuma(0:kmax,0:kmax,8,10)
complex ex,ey,ez,ex1,ey1,ez1
real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 
integer*1,pointer :: kindp(:,:) 

! set intial values; zero arrays
if (init) then
    twopi = 8.0 * atan( 1.0 )
    expsum=(0.0,0.0)
    expsuma=(0.0,0.0)
    do j=1,numvect
        sq(j,1:10)=0.0
        num(j)=0
    end do

    ! calculate number of beads of each type on each node
    nat(1:10)=0
    nata(1:10)=0
    do i=1,natms
        nat(kindp(i,1))=nat(kindp(i,1))+1
    end do

    ! collect total number of beads of each type
    call mpi_Reduce(nat, nata, 10, MPI_integer, MPI_SUM, 0 , MPI_COMM_WORLD,ierr);
    init=.false.
end if

! set basis vectors
minvec1=1.0/dlxa
minvec2=1.0/dlya
minvec3=1.0/dlza

! loop over all directions
do i=0,kmax
    do j=0,kmax
        do k=0,kmax
            
            ! check if resulting scattering vector length is less than maximum length
            len=sqrt(i*i*minvec1**2+j*j*minvec2**2+k*k*minvec3**2)*2000
            if (len>0.and. len<2000) then
                
                ! calculate number of available directions of vector projections 
                ! (if VP==0 then number of directions=1; otherwise it's=2) 
                if (i>0) then;i2st=-1;else;i2st=1;end if
                if (j>0) then;j2st=-1;else;j2st=1;end if
                if (k>0) then;k2st=-1;else;k2st=1;end if
                
                ! calculate current vector length
                rn1=minvec1*i
                rn2=minvec2*j
                rn3=minvec3*k

                ! loop over all atoms on each node considering its space position
                do atom=1,natms
                    rxb=rx(atom,1)+dlx*(ncol-1)
                    ryb=ry(atom,1)+dly*(nrow-1)
                    rzb=rz(atom,1)
                    kinda=kindp(atom,1)
                    
                    ! calculate parts of complex exponent
                    ex=cos(twopi*rn1*rxb)+(0.0,1.0)*sin(twopi*rn1*rxb)
                    ey=cos(twopi*rn2*ryb)+(0.0,1.0)*sin(twopi*rn2*ryb)
                    ez=cos(twopi*rn3*rzb)+(0.0,1.0)*sin(twopi*rn3*rzb)
                    ex1=ex
                    ey1=ey
                    ez1=ez
                    
                    ! calculate sum of exponents (at each node) considering available sqattering vector directions
                    counter=1
                    do i2=i2st,1,2
                        do j2=j2st,1,2
                            do k2=k2st,1,2

                                ! increase sum of exponents for corresponding vector lengths along y and z
                                expsum(j,k,counter,kinda)=ex1*ey1*ez1+expsum(j,k,counter,kinda) 
                                counter=counter+1
                                ez1=conjg(ez1)
                            end do
                            ey1=conjg(ey1)
                        end do
                        ex1=conjg(ex1)
                    end do
                end do
            end if
        end do
    end do
    
    ! collect sums of exponents on node 0.  do this only when scattering vector length along x is changing (to minimize number of communications)
    call mpi_Reduce(expsum, expsuma, (kmax+1)*(kmax+1)*8*10, MPI_complex, MPI_SUM, 0 , MPI_COMM_WORLD,ierr);
    
    ! add calculated sums of exponents to static structure factor (for current vector length along x and all vector lengths along y and z)
    expsum=0
    if (rank==0) then
        do j=0,kmax
            do k=0,kmax
                len=sqrt(i*i*minvec1**2+j*j*minvec2**2+k*k*minvec3**2)*2000
                if (len>0.and. len<2000) then
                    counter=1
                    if (i>0) counter=counter*2
                    if (j>0) counter=counter*2
                    if (k>0) counter=counter*2
                    do i4=1,counter
                        num(len)=num(len)+1
                        do l=1,10
                            if (nata(l)>0) sq(len,l)=sq(len,l)+conjg(expsuma(j,k,i4,l))*expsuma(j,k,i4,l)/nata(l)
                        end do
                    end do
                end if   
            end do
        end do
    end if   

! go to the next vector length along x       
end do

! output results      
if (rank==0) then
    open(10, file = 'sq.dat', status = 'unknown')
    do j=1,numvect
        if (num(j)>0) then 
            write(10,'(11f16.6,i10)')real(j)/2000,sq(j,1:10)/num(j),num(j)
        end if
    end do
    close(10)
end if

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
common /rst/ nbond,nangles
common /typel/ typelist,ntype
common /ptypeall/ kindpall
common /top/ cln,nrow,ncol 
common /angl/ ca,angles 
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
integer*4 nbond,kb,nangles
integer*1 typelist(10),ntype

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 
real*4,pointer :: vx(:),vy(:),vz(:) 
integer*1,pointer :: kindp(:,:) !
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:) 
integer*4,pointer :: cn(:,:) 
integer*1,pointer :: kindpall(:) 
integer*4,pointer :: bond(:) 
integer*4,pointer :: ca(:,:) 
integer*4,pointer :: angles(:) 
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
sz=max0(int(dlz*dly*dlx*rho*3),int(3*(2*dlz*rho*dlx*shwi+2*dlz*rho*dly*shwi+4*dlz*rho*shwi*shwi)))
ssz=max0(int(shwi*dly*dlz*11*rho),int(shwi*dlx*dlz*11*rho))
fsz=max0(int(shwi*dly*dlz*8*3*rho),int(shwi*dlx*dlz*8*3*rho))

! allocate arrays
allocate (rx(sz,3),ry(sz,3),rz(sz,3),vx(sz),vy(sz),vz(sz))
allocate (kindp(sz,3))
allocate (corr(sz,3),bcorr(npall),bcorrt(npall))
allocate (cn(npall,0:0),ca(npall,0:0))
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
ca=0
ntype=0
typelist=0

! read information about beads
do i = 1, npall
    read(1,'(i8,2i4,3f14.6)', err=10, end=10)corrb,valb,kindpb,rxb,ryb,rzb

    ! move current bead inside the simulation box
    if (rxb>dlxa) rxb=rxb-int(rxb/dlxa)*dlxa
    if (rxb<0) rxb=rxb+dlxa*(abs(int(rxb/dlxa))+1)
    if (ryb>dlya) ryb=ryb-int(ryb/dlya)*dlya
    if (ryb<0) ryb=ryb+dlya*(abs(int(ryb/dlya))+1)
    if (rzb>dlza) rzb=rzb-int(rzb/dlza)*dlza
    if (rzb<0) rzb=rzb+dlza*(abs(int(rzb/dlza))+1)
    kindpall(corrb)=kindpb
    val(corrb)=valb
  
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

! read angles information
read(1,'(a9,i9)', err=10, end=10) dummy,nangles
allocate (angles(3*nangles))
kb = 1
do i = 1, nangles
    read(1,*, err=10, end=10) i1,i2,i3
    angles(kb)=i1
    angles(kb+1)=i2
    angles(kb+2)=i3
    kb=kb+3
    ca(i1,0) = ca(i1,0) + 1
    ca(i2,0) = ca(i2,0) + 1
    ca(i3,0) = ca(i3,0) + 1
end do

! define maximum number of angles per bead
kb=maxval(ca(:,0))
nullify (ca)

! allocate angles matrix array using maximum number of bonds per bead
allocate (ca(npall,0:kb))
ca=0
kb=1

! generate angles matrix array
do i = 1, nangles
    i1=angles(kb)
    i2=angles(kb+1)
    i3=angles(kb+2)
    kb=kb+3
    ca(i1,0) = ca(i1,0) + 1
    ca(i1,ca(i1,0)) = i
    ca(i2,0) = ca(i2,0) + 1
    ca(i2,ca(i2,0)) = i
    ca(i3,0) = ca(i3,0) + 1
    ca(i3,ca(i3,0)) = i
end do
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
common /rst/ nbond,nangles
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
integer*4 kpbc,nbond,nangles
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

! allocate send/recieve arrays for bonds and angles calculation
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
common /rst/ nbond,nangles
common /bnd/ cn 
common /correl/ corr,bcorr,bcorrt
common /top/ cln,nrow,ncol 
common /angl/ ca,angles 
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
integer*4,pointer :: ca(:,:) 
integer*4,pointer :: angles(:)
integer*1,pointer :: val(:)
integer*1,pointer :: kindpall(:)

! set initial values
if (init) then
    init=.false.
    nre=0
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
    do i = 1, nangles
        write(11,'(3i8)') angles(3*i-2),angles(3*i-1),angles(3*i)
    end do
    close(11, status = 'keep')
end if

end

subroutine forcesangle()
! #################################################################
! #                                                               #
! #    subroutine 19:                                             #
! #    calculate angles forces                                    #
! #                                                               #
! #################################################################
implicit none
include "mpif.h"
save

common /rxyz/ rx,ry,rz
common /fxyz/ fx,fy,fz
common /na/ natms,natmsfan
common /cella/ dlxa,dlya,dlza
common /correl/ corr,bcorr,bcorrt
common /angl/ ca,angles
common /mpidpd/ rank, nproc
common /sizes/ sz
common /cell/ dlx,dly,dlz,rcut,kpbc
common /top/ cln,nrow,ncol 
common /prs/ press 
common /anglef/ k_eqa,k_angle
common /ptypeall/ kindpall 
common /basndrcfa/ sndb,rcvb

real*4 dlx,dly,dlz,rcut
integer*4 natms,natmsfan,mbnum,allnum
integer nam,ierr
real*4 dlxa,dlya,dlza
integer*4 i,i1,i2,j,na,k,i2b,shft,j1,j2,j3,j1a,j2a,j3a
integer*1 arrtype1,arrtype2,arrtype3,flag
integer*4 rank, nproc,cln,nrow,ncol
real*4 rxij,ryij,rzij,rijsq
real*4 w(9),g(9)
real*4 press(3)
integer*4 sz
integer*4 arrsize
integer*4 kpbc
logical init /.true./

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:)
integer*4,pointer :: mb(:) 
real*4,pointer :: fx(:),fy(:),fz(:) 
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:) 
integer*4,pointer :: cn(:,:) 
integer,pointer :: rcount(:),displ(:)
real*4,pointer :: sndb(:), rcvb(:) 
integer*4,pointer :: ca(:,:) 
integer*4,pointer :: angles(:)
real*4,pointer :: k_eqa(:,:,:),k_angle(:,:,:)
integer*1,pointer :: kindpall(:)

! allocate arrays for communications
if (init) then
    init=.false.
    allocate( mb(sz*4),rcount(nproc),displ(nproc))
    arrsize=size(sndb)
end if

! zero arrays
flag=0
mbnum=0
nam=0

! loop over all beads on each node 
do i=1,natms
    i1=corr(i,1)
    na=ca(i1,0)

    ! loop over all angles of the bead
    do k=1,na
        i2=ca(i1,k) !angle number
        j1a=angles(i2*3-2)
        j2a=angles(i2*3-1)
        j3a=angles(i2*3)
        arrtype1=bcorrt(j1a)
        arrtype2=bcorrt(j2a)
        arrtype3=bcorrt(j3a)

        ! check if there is information about all beads in the angle
        if (arrtype1*arrtype2*arrtype3/=0) then
            j1=bcorr(j1a)
            j2=bcorr(j2a)
            j3=bcorr(j3a)
            w(1) = rx(j1,arrtype1)
            w(2) = ry(j1,arrtype1)
            w(3) = rz(j1,arrtype1)
            w(4) = rx(j2,arrtype2)
            w(5) = ry(j2,arrtype2)
            w(6) = rz(j2,arrtype2)
            w(7) = rx(j3,arrtype3)
            w(8) = ry(j3,arrtype3)
            w(9) = rz(j3,arrtype3)

            ! calculate forces and internal coordinates
            call vilangle ( k_angle(kindpall(j1a),kindpall(j2a),kindpall(j3a)),k_eqa(kindpall(j1a),kindpall(j2a),kindpall(j3a)) , w, g )
            
            ! add forces and calculate virial for pressure calculation
            if(i==j1.and.arrtype1==1) then
                fx(i) = fx(i) + g(1)
                fy(i) = fy(i) + g(2)
                fz(i) = fz(i) + g(3)
            elseif (i==j2.and.arrtype2==1) then
                fx(i) = fx(i) + g(4)
                fy(i) = fy(i) + g(5)
                fz(i) = fz(i) + g(6)
                press(1)=press(1)+g(4)*w(1)+g(7)*(w(4)+w(1))
                press(2)=press(2)+g(5)*w(2)+g(8)*(w(5)+w(2))
                press(3)=press(3)+g(6)*w(3)+g(9)*(w(6)+w(3))
            elseif (i==j3.and.arrtype3==1) then
                fx(i) = fx(i) + g(7)
                fy(i) = fy(i) + g(8)
                fz(i) = fz(i) + g(9)
            end if

        ! if there is no information about some angle bead - put it to the list of lost connections
        else
            mb(mbnum+1)=i
            mb(mbnum+2)=j1a
            mb(mbnum+3)=j2a
            mb(mbnum+4)=j3a

            ! find out which bead is lost
            if (arrtype1==0) then
                nam=nam+1
                if (nam>arrsize) call error(190)
                sndb(nam)=real(j1a)
            end if
            if (arrtype2==0) then
                nam=nam+1
                if (nam>arrsize) call error(190)
                sndb(nam)=real(j2a)
            end if
            if (arrtype3==0) then
                nam=nam+1
                if (nam>arrsize) call error(190)
                sndb(nam)=real(j3a)
            end if
            mbnum=mbnum+4
        end if
    end do
end do

! gather number of lost beads
call MPI_Allgather (nam, 1, MPI_integer, rcount, 1, MPI_integer, MPI_Comm_world,ierr )
displ(1)=0
do i=2,nproc
    displ(i)=displ(i-1)+rcount(i-1)
end do
allnum=displ(nproc)+rcount(nproc)

! gather list of lost beads
call MPI_Allgatherv(sndb, nam,MPI_real4, rcvb, rcount, displ, MPI_real4, MPI_Comm_world,ierr)
nam=0

! look for information about lost beads on each node
do k=1,allnum
    if (bcorrt(int(rcvb(k)))==1) then
        i=bcorr(int(rcvb(k)))
        nam=nam+4
        if (nam>arrsize) call error(190)
        sndb(nam-3)=k
        sndb(nam-2)=rx(i,1)+(ncol-1)*dlx
        sndb(nam-1)=ry(i,1)+(nrow-1)*dly
        sndb(nam)=rz(i,1)
    end if
end do

! place information found this way to the array in the order of receiving
rcvb(1:allnum*3)=0
do k=1,nam,4
    i=sndb(k)*3-2
    rcvb(i)=sndb(k+1)
    rcvb(i+1)=sndb(k+2)
    rcvb(i+2)=sndb(k+3)
end do

! gather linformation about found angles beads
call mpi_allreduce(mpi_in_place,rcvb,allnum*3,mpi_real,mpi_sum,mpi_comm_world,ierr)

! use collected information to calculate angle forces
i1=displ(rank+1)*3+1
do j=1,mbnum,4
    i=mb(j)
    j1a=mb(j+1)
    j2a=mb(j+2)
    j3a=mb(j+3)
    arrtype1=bcorrt(j1a)
    arrtype2=bcorrt(j2a)
    arrtype3=bcorrt(j3a)

    ! get beads coordinates
    if (arrtype1==0) then
        w(1)=rcvb(i1)-(ncol-1)*dlx
        w(2)=rcvb(i1+1)-(nrow-1)*dly
        w(3)=rcvb(i1+2)
        i1=i1+3
    else
        j1=bcorr(j1a)
        w(1) = rx(j1,arrtype1)
        w(2) = ry(j1,arrtype1)
        w(3) = rz(j1,arrtype1)
    end if
    if (arrtype2==0) then
        w(4)=rcvb(i1)-(ncol-1)*dlx
        w(5)=rcvb(i1+1)-(nrow-1)*dly
        w(6)=rcvb(i1+2)
        i1=i1+3
    else
        j2=bcorr(j2a)
        w(4) = rx(j2,arrtype2)
        w(5) = ry(j2,arrtype2)
        w(6) = rz(j2,arrtype2)
    end if
    if (arrtype3==0) then
        w(7)=rcvb(i1)-(ncol-1)*dlx
        w(8)=rcvb(i1+1)-(nrow-1)*dly
        w(9)=rcvb(i1+2)
        i1=i1+3
    else
        j3=bcorr(j3a)
        w(7) = rx(j3,arrtype3)
        w(8) = ry(j3,arrtype3)
        w(9) = rz(j3,arrtype3)
    end if

    ! calculate forces and internal coordinates
    call vilangle ( k_angle(kindpall(j1a),kindpall(j2a),kindpall(j3a)),k_eqa(kindpall(j1a),kindpall(j2a),kindpall(j3a)), w, g )

    ! calculate forces and virial for pressure calculation
    if(i==j1.and.arrtype1==1) then
        
        ! if 1st bead in the angle is a real bead
        fx(i) = fx(i) + g(1)
        fy(i) = fy(i) + g(2)
        fz(i) = fz(i) + g(3)
    elseif (i==j2.and.arrtype2==1) then

        ! if 2nd bead in the angle is a real bead
        fx(i) = fx(i) + g(4)
        fy(i) = fy(i) + g(5)
        fz(i) = fz(i) + g(6)
        press(1)=press(1)+g(4)*w(1)+g(7)*(w(4)+w(1))
        press(2)=press(2)+g(5)*w(2)+g(8)*(w(5)+w(2))
        press(3)=press(3)+g(6)*w(3)+g(9)*(w(6)+w(3))
    elseif (i==j3.and.arrtype3==1) then

        ! if 3rd bead in the angle is a real bead
        fx(i) = fx(i) + g(7)
        fy(i) = fy(i) + g(8)
        fz(i) = fz(i) + g(9)
    end if
end do

end

subroutine vilangle( ca, a0, x, f )
! #################################################################
! #                                                               #
! #    subroutine 20:                                             #
! #    calculate forces acting on angle beads using               #
! #    vilson s-vector 1-2-3                                      #
! #                                                               #
! #################################################################
implicit real*4 (a-h,o-z)
real*4 x(9),f(9),e1(3),e2(3)

call one(x(4),x(1),e1,b1)
call one(x(7),x(4),e2,b2)

cosa=-dotd(e1,e2)
qimj=ang(cosa)

sina=sqrt(abs(1.0-cosa*cosa))

if(sina<1.0e-5)go to 20
da=(qimj-a0)*1.7453292e-2
a1=1.0/(b1*sina)
a2=-1.0/(b2*sina)
uq=ca*da*da
f0=-2.0*ca*da
do 10 k=1,3
    s1=(e1(k)*cosa+e2(k))*a1
    s3=(e2(k)*cosa+e1(k))*a2
    s2=-s1-s3
    f(k)=f0*s1
    f(k+3)=f0*s2
    f(k+6)=f0*s3
10 continue
return
20 uq=0.0
do 30 i=1,9
30 f(i)=0.0

return
end

subroutine one(x1,x2,e,r)
! #################################################################
! #                                                               #
! #    subroutine 21:                                             #
! #    calculate unit vector and internal coordinates of beads    #
! #                                                               #
! #################################################################
implicit real*4 (a-h,o-z)
common /cella/ dlxa,dlya,dlza    
real*4 x1(3),x2(3),e(3)
real*4 dlxa,dlya,dlza

a=x2(1)-x1(1)
b=x2(2)-x1(2)
c=x2(3)-x1(3)

! apply PBC
if (a>dlxa/2) then 
    a=a-dlxa
elseif (a<-dlxa/2) then
    a=a+dlxa
end if      
if (b>dlya/2) then 
    b=b-dlya
elseif (b<-dlya/2) then
    b=b+dlya
end if    
if (c>dlza/2) then 
    c=c-dlza
elseif (c<-dlza/2) then
    c=c+dlza
end if

! calculate vector length
r=sqrt(a*a+b*b+c*c)+1.0e-5

! calculate angle internal coordinates
x2(1)=-a
x2(2)=-b
x2(3)=-c

! calcuclate unit vector
e(1)=a/r
e(2)=b/r
e(3)=c/r

end

real*4 function ang ( c )
! #################################################################
! #                                                               #
! #    function 22:                                               #
! #    find angle (in grad), using its cosine (c)                 #
! #                                                               #
! #################################################################
implicit real*4 (a-h, o-z)
save

data key/0/

if(key/=0)go to 20
pi=4.0*atan(1.e0)
d=180.0/pi
e=1.0
10 e=e/2.0
t=1.0 + e
if(t>1.0)go to 10
e=sqrt(e)
key=1

20 a=pi/2.0
if(abs(c)<=0.0)go to 30
a=atan(sqrt(abs(1.0-c*c))/c)
if(a<0.0)a=a+pi
30 ang=a*d

end

real*4 function dotd(e1,e2)
! #################################################################
! #                                                               #
! #    function 23:                                               #
! #    calculate scalar product of vectors                        #
! #                                                               #
! #################################################################
implicit real*4 (a-h,o-z)
real*4 e1(3),e2(3)

dotd=e1(1)*e2(1)+e1(2)*e2(2)+e1(3)*e2(3)

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

subroutine deform(tau,tai)
! #################################################################
! #                                                               #
! #    subroutine 25:                                             #
! #    perform affine deformation of the system                   #
! #                                                               #
! #################################################################
implicit none
save

common /rxyz/ rx,ry,rz
common /cella/ dlxa,dlya,dlza
common /na/ natms,natmsfan
common /cell/ dlx,dly,dlz,rcut,kpbc
common /bxsize/ ilx,ily,ilz

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:)

real*8 mu(3),volm
real*4 dlxa,dlya,dlza
real*4 dlx,dly,dlz,rcut
integer*4 natms,natmsfan,kpbc,ilx,ily,ilz
real*4 tai,tau
logical init /.true./

! calculate initial values
if(init) then
    volm = dlxa*dlya*dlza
    mu(1)=tai**2
    mu(2)=tai
    mu(3)=tai
    init=.false.
end if

! stretching along x
if (tau>=0) then
    dlxa=dlxa*mu(1)
    dlya=dlya/mu(2)
    dlza=dlza/mu(3)
    dlx=dlx*mu(1)
    dly=dly/mu(2)
    dlz=dlz/mu(3)
    ilx=int(dlx/rcut)
    ily=int(dly/rcut)
    ilz=int(dlz/rcut)
    rx(1:natms,1)=rx(1:natms,1)*mu(1)
    ry(1:natms,1)=ry(1:natms,1)/mu(2)
    rz(1:natms,1)=rz(1:natms,1)/mu(3)
    if (dly<rcut) call error(251)
    if (dlz<rcut) call error(252)
! compression along x
else
    dlxa=dlxa/mu(1)
    dlya=dlya*mu(2)
    dlza=dlza*mu(3)
    dlx=dlx/mu(1)
    dly=dly*mu(2)
    dlz=dlz*mu(3)
    ilx=int(dlx/rcut)
    ily=int(dly/rcut)
    ilz=int(dlz/rcut)
    rx(1:natms,1)=rx(1:natms,1)/mu(1)
    ry(1:natms,1)=ry(1:natms,1)*mu(2)
    rz(1:natms,1)=rz(1:natms,1)*mu(3)
    if (dlx<rcut) call error(253)
end if

end

subroutine pressout(flag)
! #################################################################
! #                                                               #
! #    subroutine 26:                                             #
! #    output stress-strain curve                                 #
! #                                                               #
! #################################################################
implicit none
save

logical init /.true./
common /cella/ dlxa,dlya,dlza
common /prs/ press
common /defor/ nave 

integer*4 steps3,qstp,nout,nit,i,nit2
real*4 dlxa,dlya,dlza,L0
real*4 press(3)
integer flag
integer*4 nave
real*8 volm,sd
real*8,pointer :: presx(:),presy(:),presz(:)

! allocate pressure arrays, set initial values
if (init) then
    init=.false.
    allocate(presx(nave),presy(nave),presz(nave))
    presx=0
    presy=0
    presz=0
    sd=0
    L0=dlxa
    open (1,file='Stress-Strain.dat')
    close (1,status='delete')
    nit=0
    volm = dlxa*dlya*dlza
end if

! collect pressures for averaging
if (flag==0) then
    nit=nit+1
    presx(nit)=dble(press(1))
    presy(nit)=dble(press(2))
    presz(nit)=dble(press(3))

! output stress-strain curve
else
    open (1,file='Stress-Strain.dat',position= 'append')

    ! calculate standard deviation
    if (nit>1) then
        do i=1,nit
            sd=sd+(presx(i)-sum(presx(1:nit))/dble(nit))**2
        end do
        sd=sqrt(sd/(nit-1))
    end if
    write(1,'(9f16.9,i8)'),dlxa/L0,sum(presx(1:nit))/dble(nit),sd,sum(presy(1:nit))/dble(nit),sum(presz(1:nit))/dble(nit),dlxa,dlya,dlza,dlxa*dlya*dlza/volm-1,nit
    nit=0
    close(1)
end if

end

subroutine chem()
! #################################################################
! #                                                               #
! #    subroutine 27:                                             #
! #    chemical reactions                                         #
! #                                                               #
! #################################################################
implicit none
save
include "mpif.h"

common /bxsize/ ilx,ily,ilz
common /cell/ dlx,dly,dlz,rcut,kpbc
common /rst/ nbond,nangles
common /ptype/ kindp 
common /npall/ npall 
common /bnd/ cn
common /angl/ ca,angles
common /sizes/ sz 
common /mpidpd/ rank, nproc 
common /correl/ corr,bcorr,bcorrt 
common /ptypeall/ kindpall 
common /top/ cln,nrow,ncol 
common /chemmap/ nmap,arrtypemap
common /probcb/ probc,probb 
common /vale/ val 
common /ch/ rchem 


integer*1,pointer :: kindp(:,:) 
integer*4,pointer :: cn(:,:) 
integer*4,pointer :: corr(:,:),bcorr(:)
integer*1,pointer :: bcorrt(:) 
integer*1,pointer :: kindpall(:) 
integer,pointer :: rcount(:),displ(:) 
integer*4,pointer :: reacs(:), rcvb(:) 
integer*4,pointer :: nmap(:,:,:,:,:) 
integer*1,pointer :: arrtypemap(:,:,:,:) 
integer*1,pointer :: val(:) 
real*4,pointer :: probc(:,:),probb(:,:)
integer*4,pointer :: ca(:,:) 
integer*4,pointer :: angles(:)

integer*4 rank,nproc,cln,nrow,ncol,rowd,cold,rowc,colc
integer*4 nbond,nangles
integer*4 sz
integer*4 link(sz,3),lct(0:ilx+1,0:ily+1,ilz),qua(0:ilx+1,0:ily+1,ilz)
real*4 dlx,dly,dlz,rcut
integer*4 ilx,ily,ilz,kpbc,ntry
real*4 rchem
integer*4 i,j,k,npall,iter,iat,ix,iy,iz,klast,l,ic,jc,ineig,jx,jy,jz,iatm1,iatm2,nami,namj,vali,valj
real*4 uni
integer*4 reaclist(1000),nneib,reaclistc(1000),reaclistart(1000)
real*4 distlist(1000),rij,rchem2
integer*4 arrtype
integer*4 reacssize,allnum
integer ierr
logical init /.true./
logical change

! sey some initial values
if (init) then
    rchem=rcut

    ! allocate main chem arrays
    allocate( rcount(nproc), displ(nproc),reacs(npall),rcvb(npall*2),nmap(3,27,ilx,ily,ilz),arrtypemap(27,ilx,ily,ilz))
    
    ! create map of subcells
    call maps()
    
    rchem2=rchem**2

    ! divide nodes into groups
    cold=min0(3,cln)
    rowd=min0(3,nproc/cln)
    if (dlx>rchem*2) cold=min0(2,cln)
    if (dly>rchem*2) rowd=min0(2,nproc/cln)

    ! number of iterations at every choice of nodes group
    ntry=real(npall*cold*rowd)/real(200*nproc)
    init=.false.
end if

! create linked lists for all beads (for real and ghost)
call links(link,qua,lct)

! change groups of nodes
do iter=1,200
    reacssize=0 

    ! choose nodes group at random
    if (rank==0) then
        rowc=int(uni()*rowd)
        colc=int(uni()*cold)
    end if
    call MPI_Bcast ( rowc, 1, MPI_integer, 0,mpi_comm_world,ierr )
    call MPI_Bcast ( colc, 1, MPI_integer, 0,mpi_comm_world,ierr ) 
    
    ! if the node belongs to the chosen group        
    if (mod(nrow-1,rowd)==rowc.and.mod(ncol-1,cold)==colc) then
    
        ! choose ntry beads at random 
        do iat=1,ntry
            
            ! choose a bead at random
            ix=int(uni()*ilx)+1
            iy=int(uni()*ily)+1
            iz=int(uni()*ilz)+1
            klast=qua(ix,iy,iz)
            if (klast==0) cycle
            k=int(uni()*klast)
            i=lct(ix,iy,iz)
            do l=1,k
                i=link(i,1)
            end do
            ic=corr(i,1)
            vali=val(ic)
            
            ! check if its valency is not equal to 0    
            if (vali>0) then
                nneib=0  

                ! loop over all neighbouring subcells
                do ineig=1,27
                    jx=nmap(1,ineig,ix,iy,iz)
                    jy=nmap(2,ineig,ix,iy,iz)
                    jz=nmap(3,ineig,ix,iy,iz)
                    arrtype=arrtypemap(ineig,ix,iy,iz)

                    ! choose the first bead in the selected subcell
                    j=lct(jx,jy,jz)

                    ! loop over all beads in the selected subcell
                    do while (j/=0)
                        if ((i==j).and.(arrtype==1)) then
                            j=link(j,arrtype)
                            cycle
                        end if
                        
                        ! check if the valency of j bead is not equal to 0
                        jc=corr(j,arrtype)
                        valj=val(jc)
                        if (valj<1) then
                            j=link(j,arrtype)
                            cycle
                        end if
                        
                        ! check distance between beads
                        call distij(i,j,arrtype,rij)
                        if (rij>=rchem2) then
                            j=link(j,arrtype)
                            cycle
                        end if
                        
                        ! place j bead to the list of possible connections
                        nneib=nneib+1
                        reaclistc(nneib)=jc
                        reaclist(nneib)=j
                        reaclistart(nneib)=arrtype
                        distlist(nneib)=rij
                        j=link(j,arrtype)
                    end do
                end do
                
                ! sort possible connections list ascending
                call sort(nneib,reaclist,reaclistc,reaclistart,distlist)
                
                ! loop over list of possible connections
                do ineig=1,nneib
                    if (vali<1) exit
                    j=reaclist(ineig)
                    arrtype=reaclistart(ineig)

                    ! probability test
                    if (abs(corr(i,1)-corr(j,arrtype))>1 .and. uni()<probc(kindp(i,1),kindp(j,arrtype))) then

                        ! form a bond
                        jc=reaclistc(ineig)
                        cn(ic,0) = cn(ic,0)+1
                        cn(ic,cn(ic,0)) = jc
                        cn(jc,0) = cn(jc,0)+1
                        cn(jc,cn(jc,0)) = ic
                        val(ic)=val(ic)-1
                        val(jc)=val(jc)-1
                        vali=vali-1
                        reacs(reacssize+1)=ic
                        reacs(reacssize+2)=jc
                        reacssize=reacssize+2
                        nbond=nbond+1
                    end if
                end do
            end if
        end do
    end if

    ! gather number of bonds created on all nodes
    call MPI_Allgather (reacssize, 1, MPI_integer, rcount, 1, MPI_integer, mpi_comm_world,ierr )
    displ(1)=0
    do i=2,nproc
        displ(i)=displ(i-1)+rcount(i-1)
    end do
    allnum=displ(nproc)+rcount(nproc)    

    ! receive information about all created bonds
    call MPI_Allgatherv(reacs, reacssize,MPI_integer, rcvb, rcount, displ, MPI_integer, mpi_comm_world,ierr)

    ! on each node create bonds (add them to the connection matrix) obtained from other nodes
    do i=0,displ(rank+1)-1,2
        ic=rcvb(i+1)
        jc=rcvb(i+2)
        cn(ic,0)=cn(ic,0)+1
        cn(jc,0)=cn(jc,0)+1
        cn(ic,cn(ic,0))=jc
        cn(jc,cn(jc,0))=ic
        val(ic)=val(ic)-1
        val(jc)=val(jc)-1
        nbond=nbond+1
    end do
    if (rank+1<nproc) then
        do i=displ(rank+2),allnum-1,2
            ic=rcvb(i+1)
            jc=rcvb(i+2)
            cn(ic,0)=cn(ic,0)+1
            cn(jc,0)=cn(jc,0)+1
            cn(ic,cn(ic,0))=jc
            cn(jc,cn(jc,0))=ic
            val(ic)=val(ic)-1
            val(jc)=val(jc)-1
            nbond=nbond+1
        end do
    end if    
end do

! bond breaking part
! divide all beads among all nodes
iatm1 = 1+(rank*npall)/nproc
iatm2 =  ((rank+1)*npall)/nproc
reacssize=0

! loop over beads
do iat=iatm1,iatm2
    do i=1,cn(iat,0)
        if(iat<cn(iat,i)) then
            
            ! angle bonds can't be broken
            change=.false.
            do j=1,ca(iat,0)
                if(angles(ca(iat,j)*3-2)==cn(iat,i).or.angles(ca(iat,j)*3-1)==cn(iat,i).or.angles(ca(iat,j)*3)==cn(iat,i)) change=.true.
            end do

            ! probability test 
            if (abs(cn(iat,i)-iat)>1 .and. uni()<probb(kindpall(iat),kindpall(cn(iat,i)))) then
                if (change) cycle
                reacs(reacssize+1)=iat
                reacs(reacssize+2)=i
                reacssize=reacssize+2
            end if
        end if  
    end do
end do

! gather number of destroyed bonds
call MPI_Allgather (reacssize, 1, MPI_integer, rcount, 1, MPI_integer, mpi_comm_world,ierr )
displ(1)=0
do i=2,nproc
    displ(i)=displ(i-1)+rcount(i-1)
end do
allnum=displ(nproc)+rcount(nproc)    

! gather all information about broken bonds
call MPI_Allgatherv(reacs, reacssize,MPI_integer, rcvb, rcount, displ, MPI_integer, mpi_comm_world,ierr)

! remove all broken bonds from the connection matrix
do i=0,allnum-1,2
    nbond=nbond-1
    ic=rcvb(i+1)
    k=rcvb(i+2)
    jc=cn(ic,k)
    val(ic)=val(ic)+1
    val(jc)=val(jc)+1
    cn(ic,k)=cn(ic,cn(ic,0))
    cn(ic,0)=cn(ic,0)-1
    do j = 1, cn(jc,0)
        if (cn(jc,j)==ic)then
            cn(jc,j) = cn(jc,cn(jc,0))
            cn(jc,0) = cn(jc,0) - 1
            exit
        end if
    end do
end do

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
ryb=ryb-dlya*nint(rzb/dlya)
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

subroutine links(link,qua,lct)
! #################################################################
! #                                                               #
! #    subroutine 30:                                             #
! #    create linked lists for all beads                          #
! #                                                               #
! #################################################################  
implicit none
save

common /cell/ dlx,dly,dlz,rcut,kpbc
common /ch/ rchem 
common /na/ natms,natmsfan
common /nab/ natmsfanbond
common /sizes/ sz
common /bxsize/ ilx,ily,ilz 
common /rxyz/ rx,ry,rz 

real*4,pointer :: rx(:,:),ry(:,:),rz(:,:) 

real*4 dlx,dly,dlz,rcut
integer*4 sz
integer*4 ilx,ily,ilz,kpbc
real*4 rchem
integer*4 link(sz,3),lct(0:ilx+1,0:ily+1,ilz),qua(0:ilx+1,0:ily+1,ilz)
integer*4 natms,natmsfan,natmsfanbond
integer*4 i,j,ix,iy,iz

! zero arrays
link=0
qua=0
lct=0

! create linked list for real beads
do i=1,natms
    ix=min(int(rx(i,1)/rchem),ilx-1)+1
    iy=min(int(ry(i,1)/rchem),ily-1)+1
    iz=min(int(rz(i,1)/rchem),ilz-1)+1
    j=lct(ix,iy,iz)
    lct(ix,iy,iz)=i
    link(i,1)=j
    qua(ix,iy,iz)=qua(ix,iy,iz)+1
end do

! create linked list for ghost beads type 1(which was used to calculate forces)
do i=1,natmsfan
    if (rx(i,2)>=(dlx+rchem)) then
        cycle
    elseif (rx(i,2)<dlx) then
        ix=min(int(rx(i,2)/rchem),ilx-1)+1
    elseif (rx(i,2)<(dlx+rchem)) then
        ix=ilx+1
    end if
    if (ry(i,2)>=(dly+rchem).or.ry(i,2)<=(-rchem)) then
        cycle
    elseif (ry(i,2)<0) then
        iy=0  
    elseif (ry(i,2)<dly) then
        iy=min(int(ry(i,2)/rchem),ily-1)+1
    elseif (ry(i,2)<(dly+rchem)) then
        iy=ily+1
    end if
    iz=min(int(rz(i,2)/rchem),ilz-1)+1
    j=lct(ix,iy,iz)
    lct(ix,iy,iz)=i
    link(i,2)=j
    qua(ix,iy,iz)=qua(ix,iy,iz)+1
end do

! create linked list for ghost beads type 2(which was not used to calculate forces)
do i=1,natmsfanbond
    if (rx(i,3)<=-(rchem)) then
        cycle
    elseif (rx(i,3)<=0) then
        ix=0
    elseif (rx(i,3)<=dlx) then
        ix=min(int(rx(i,3)/rchem),ilx-1)+1
    end if
    if (ry(i,3)>=(dly+rchem).or.ry(i,3)<=(-rchem)) then
        cycle
    elseif (ry(i,3)<=0) then
        iy=0  
    elseif (ry(i,3)<dly) then
        iy=min(int(ry(i,3)/rchem),ily-1)+1
    elseif (ry(i,3)<(dly+rchem)) then
        iy=ily+1
    end if
    iz=min(int(rz(i,3)/rchem),ilz-1)+1
    j=lct(ix,iy,iz)
    lct(ix,iy,iz)=i
    link(i,3)=j
    qua(ix,iy,iz)=qua(ix,iy,iz)+1
end do    

end

subroutine maps()
! #################################################################
! #                                                               #
! #    subroutine 31:                                             #
! #    create full map of cubcells neighbours                     #
! #    (all 27 unlike fmaps)                                      #
! #    for chem                                                   #
! #                                                               #
! #################################################################  
implicit none

common /bxsize/ ilx,ily,ilz
common /chemmap/ nmap,arrtypemap

integer*1,pointer :: arrtypemap(:,:,:,:)
integer*4,pointer :: nmap(:,:,:,:,:)

integer*4 ilx,ily,ilz
integer*4 ix,iy,iz,jx,jy,jz,nei
integer*4 nix,niy,niz
dimension nix(27),niy(27),niz(27)
data niz/0,1,0,0,-1,1, 0,-1,1,0,-1, 1,-1,1,-1, 0, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1/
data niy/0,0,1,0, 1,1,-1, 0,0,1,-1,-1, 1,1, 0,-1, 0,-1,-1, 1, 0, 0,-1, 1, 1,-1,-1/
data nix/0,0,0,1, 0,0, 1, 1,1,1, 1, 1, 1,1, 0, 0,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1/

! loop over all subcells
do ix=1,ilx
    do iy=1,ily
        do iz=1,ilz

            ! loop over all possible neighbouring subcells
            do nei=1,27
                jx=ix+nix(nei)
                jy=iy+niy(nei)
                jz=iz+niz(nei)
                
                ! apply PBC
                if(jz>ilz)then
                    jz=1
                elseif(jz<1)then
                    jz=ilz
                end if

                ! update subcells map arrays
                nmap(1,nei,ix,iy,iz)=jx
                nmap(2,nei,ix,iy,iz)=jy
                nmap(3,nei,ix,iy,iz)=jz

                ! create map of ghost subcells type
                if (jx>ilx) then
                    arrtypemap(nei,ix,iy,iz)=2
                elseif (jx<1) then
                    arrtypemap(nei,ix,iy,iz)=3
                elseif (jy>ily) then
                    arrtypemap(nei,ix,iy,iz)=2
                elseif (jy<1) then
                    arrtypemap(nei,ix,iy,iz)=3
                else
                    arrtypemap(nei,ix,iy,iz)=1              
                end if
            end do
        end do
    end do
end do

end

subroutine rconf(steps1,steps2,steps3,qstp,rst,rststp,vtk,vtkstp,sqf,sqfstp,chemstp,stepsdef,ndef,tau,tai,def_res)
! #################################################################
! #                                                               #
! #    subroutine 32:                                             #
! #    read input script                                          #
! #                                                               #
! #################################################################  
implicit none

common /cella/ dlxa,dlya,dlza 
common /den/ rho 
common /data/ alpha,sigma,gamma
common /time/ dt
common /shwi/ shwi
common /probcb/ probc,probb
common /mpidpd/ rank, nproc 
common /bondf/ k_eq,k_bond 
common /anglef/ k_eqa,k_angle
common /bmult/ mult
common /defor/ nave 
common /lambd/ lambda

character*16,pointer :: fnames(:)
real*4,pointer :: conc(:)
real*4,pointer :: alpha(:,:),probc(:,:),probb(:,:)
integer*1,pointer :: anglestempl(:,:)
integer*1,pointer :: vallist(:),place(:)
character*2,pointer :: bead_type(:)
real*4,pointer :: k_eq(:,:),k_bond(:,:)
real*4,pointer :: k_eqa(:,:,:),k_angle(:,:,:)

character*256 str
character*2 el

integer*4 i,j,k
integer*1 rtype,snum,spnum,s1,s2,atempl
real*4 dlxa,dlya,dlza,dt,sigma,gamma
integer*4 chemstp, steps1,steps2,steps3,stepsdef
real*4 alph,prob,probr
real*4 rho,shwi
integer*4 rst,vtk,sqf,qstp,vtkstp,sqfstp,rststp,nave,ndef
real*4 tai,tau
integer*4 rank,nproc,mult
logical exists
real*4 graftdens
integer*4 ntlength
real*4 lambda
integer*1 def_res
real*4 sphrad

! set default values
rst=0
vtk=0
sqf=0
chemstp=200
sigma=3.0
shwi=1.0
vtkstp=10000
rststp=50000
steps2=0
steps1=0
steps3=0
stepsdef=0
spnum=0
mult=1
dt=0.04
qstp=1000
ntlength=40
sphrad=2.0
graftdens=0.0
def_res=0

! check if input script exists
inquire( file = 'dpdconf.dat', exist = exists )
if ( .not. exists ) goto 14

open(1, file = 'dpdconf.dat')
read (1,'(a256)') str
do i=1,256
    if (str(i:i)/=' '.and.str(i:i)/=char(09)) exit
end do
rtype=0

! read type of simulation
if (str(i:i+len('run_type')-1)=='run_type') then
    read (str(i+len('run_type'):256),*,err=15,end=15) el
    if (el(1:2)=='co') then
        rtype=1
    elseif (el(1:2)=='re') then
        rtype=2
    else
        call error(320)
    end if
end if

atempl=0

! allocate all arrays
allocate(k_eq(10,10),k_bond(10,10),alpha(10,10),probc(10,10),anglestempl(1000,3),vallist(10),probb(10,10),bead_type(10),k_eqa(10,10,10),k_angle(10,10,10))
vallist=0
probc=0
probb=0
alpha=25.0
k_eq=0.0
k_eqa=180.0
k_bond=4.0
k_angle=0.0
snum=0
bead_type='  '
lambda=0.65

! read information from script file
do
    if (rtype==0) call error(321)

    read (1,'(a256)') str
    do i=1,256
        if (str(i:i)/=' '.and.str(i:i)/=char(09)) exit
    end do
    if (i==257) cycle

    if (str(i:i+2)=='end') exit

    if (str(i:i+9)=='input_file') then
        j=i+10
        if (snum==0.and.rtype==1) call error(322)
        do
            do i=j,256
                if (str(i:i)/=' '.and.str(i:i)/=char(09)) exit
            end do
        
            if (i==257) exit
            spnum=spnum+1
            k=1
            do j=i,256
                if (str(j:j)/=' '.and.str(j:j)/=char(09))then
                    fnames(spnum)(k:k)=str(j:j)
                    k=k+1
                else
                    exit
                end if
            end do
        end do
    end if

    if (i==257) cycle

    if (str(i:i+len('species_number')-1)=='species_number') then
        read (str(i+len('species_number'):256),*,err=4,end=4) snum
        allocate (fnames(snum),conc(snum),place(snum))
        place=1
        fnames=' '
        conc=0
    end if
    
    if (str(i:i+13)=='concentrations') then
        if (snum==0) goto 1
        read (str(i+14:256),*,err=2,end=2) conc(1:snum)
        if (sum(conc(1:snum))/=100.and.rtype==1) goto 3
    end if

    if (str(i:i+len('place')-1)=='place') then
        if (snum==0) goto 1
        read (str(i+len('place'):256),*,err=21,end=21) s1,s2
        place(s1)=s2
    end if

    if (str(i:i+7)=='box_size') then
        read (str(i+8:256),*,err=5,end=5) dlxa,dlya,dlza
    end if

    if (str(i:i+7)=='timestep') then
        read (str(i+8:256),*,err=7,end=7) dt
    end if

    if (str(i:i+6)=='steps_1') then
        read (str(i+7:256),*,err=8,end=8) steps1
    end if

    if (str(i:i+6)=='steps_2') then
        read (str(i+7:256),*,err=9,end=9) steps2
    end if

    if (str(i:i+6)=='steps_3') then
        read (str(i+7:256),*,err=10,end=10) steps3
    end if

    if (str(i:i+7)=='steps_ch') then
        read (str(i+8:256),*,err=11,end=11) chemstp
    end if

    if (str(i:i+9)=='pair_coeff') then
        read (str(i+10:256),*,err=12,end=12) s1,s2,alph

        ! apply wildcards
        if (s1==0.and.s2/=0) then
            do s1=1,10
                alpha(s1,s2)=alph
                alpha(s2,s1)=alph
            end do
        elseif (s1/=0.and.s2==0) then
            do s2=1,10
                alpha(s1,s2)=alph
                alpha(s2,s1)=alph
            end do
        elseif (s1==0.and.s2==0) then
            do s2=1,10
                do s1=1,10
                    alpha(s1,s2)=alph
                    alpha(s2,s1)=alph
                end do
            end do
        else
            alpha(s1,s2)=alph
            alpha(s2,s1)=alph
        end if
    end if

    if (str(i:i+9)=='chem_coeff') then
        read (str(i+10:256),*,err=13,end=13) s1,s2,prob,probr
        probc(s1,s2)=prob
        probc(s2,s1)=prob
        probb(s1,s2)=probr
        probb(s2,s1)=probr
    end if

    if (str(i:i+len('dpd_density')-1)=='dpd_density') then
        read (str(i+len('dpd_density'):256),*,err=15,end=15) rho
    end if

    if (str(i:i+len('output_freq')-1)=='output_freq') then
        read (str(i+len('output_freq'):256),*,err=16,end=16) qstp
    end if

    if (str(i:i+len('angle_template')-1)=='angle_template') then
        atempl=atempl+1
        read (str(i+len('angle_template'):256),*,err=17,end=17) anglestempl(atempl,1),anglestempl(atempl,2),anglestempl(atempl,3),alph,prob
        k_eqa(anglestempl(atempl,1),anglestempl(atempl,2),anglestempl(atempl,3))=alph
        k_eqa(anglestempl(atempl,3),anglestempl(atempl,2),anglestempl(atempl,1))=alph
        k_angle(anglestempl(atempl,3),anglestempl(atempl,2),anglestempl(atempl,1))=prob
        k_angle(anglestempl(atempl,1),anglestempl(atempl,2),anglestempl(atempl,3))=prob
    end if

    if (str(i:i+len('bond_template')-1)=='bond_template') then
        read (str(i+len('bond_template'):256),*,err=18,end=18) s1,s2,alph,prob

        ! apply wildcards
        if (s1==0.and.s2/=0) then
            do s1=1,10
                k_eq(s1,s2)=alph
                k_eq(s2,s1)=alph
                k_bond(s1,s2)=prob
                k_bond(s2,s1)=prob
            end do
        elseif (s1/=0.and.s2==0) then
            do s2=1,10
                k_eq(s1,s2)=alph
                k_eq(s2,s1)=alph
                k_bond(s1,s2)=prob
                k_bond(s2,s1)=prob
            end do
        elseif (s1==0.and.s2==0) then
            do s2=1,10
                do s1=1,10
                    k_eq(s1,s2)=alph
                    k_eq(s2,s1)=alph
                    k_bond(s1,s2)=prob
                    k_bond(s2,s1)=prob
                end do
            end do
        else
            k_eq(s1,s2)=alph
            k_eq(s2,s1)=alph
            k_bond(s1,s2)=prob
            k_bond(s2,s1)=prob
        end if
    end if

    if (str(i:i+len('shwi')-1)=='shwi') then
        read (str(i+len('shwi'):256),*,err=26,end=26) shwi
    end if

    if (str(i:i+len('valency')-1)=='valency') then
        read (str(i+len('valency'):256),*,err=19,end=19) s1,s2
        vallist(s1)=s2
    end if

    if (str(i:i+len('deform')-1)=='deform') then
        read (str(i+len('deform'):256),*,err=20,end=20) stepsdef,ndef,tau,tai,nave
    end if

    

    if (str(i:i+len('sqf')-1)=='sqf') then
        read (str(i+len('sqf'):256),*,err=22,end=22) sqfstp
        sqf=1
    end if

    if (str(i:i+len('restart')-1)=='restart') then
        read (str(i+len('restart'):256),*,err=23,end=23) rststp
        rst=1
    end if

    if (str(i:i+len('vtk')-1)=='vtk') then
        read (str(i+len('vtk'):256),*,err=24,end=24) vtkstp
        vtk=1
    end if

    if (str(i:i+len('bead_type')-1)=='bead_type') then
        read (str(i+len('bead_type'):256),*,err=25,end=25) s1,el
        if(s1>10) goto 25
        bead_type(s1)=el
    end if

    if (str(i:i+len('nt_shape')-1)=='nt_shape') then
        read (str(i+len('nt_shape'):256),*,err=27,end=27) ntlength,graftdens
    end if

    if (str(i:i+len('lambda')-1)=='lambda') then
        read (str(i+len('lambda'):256),*,err=28,end=28) lambda
    end if

    if (str(i:i+len('def_res')-1)=='def_res') then
        read (str(i+len('def_res'):256),*,err=29,end=29) def_res
    end if

    if (str(i:i+len('sph_shape')-1)=='sph_shape') then
        read (str(i+len('sph_shape'):256),*,err=30,end=30) sphrad,graftdens
    end if
end do

gamma=0.50*sigma*sigma

! construct initial system
if (rtype==1.and.rank==0) call construct(fnames,conc,place,snum,anglestempl,atempl,vallist,bead_type,ntlength,sphrad,graftdens)
call barrier

! delete temporary arrays
nullify(fnames,conc,anglestempl,vallist,place)

return
1 call error(323)    !stop 'error:spnum=0'
2 call error(324)    !    stop 'error:conc number'
3 call error(325)    !    stop 'error:conc/=100'
4 call error(326)    !    stop 'error:n_species error'
5 call error(327)    !    stop 'error:box size error'
6 call error(328)    !    stop 'error:prname error'
7 call error(329)    !    stop 'error:dt error'
8 call error(3210)    !    stop 'error:t1 error'
9 call error(3211)    !    stop 'error:t1 error'
10 call error(3212)    !   stop 'error:t1 error'
11 call error(3213)    !   stop 'error:tch error'
12 call error(3214)    !   stop 'error:alpha error'
13 call error(3215)    !   stop 'error:probc error'
14 call error(3216)    !   stop 'there is no dpdconf.dat file! '
15 call error(3217)    !   stop 'error:density error!'
16 call error(3218)    !   stop 'output_freq error'
17 call error(3219)    !   stop 'angletempl error'
18 call error(3220)    !   stop 'btempl error'
19 call error(3221)    !   stop 'val error'
20 call error(3222)    !   stop 'deform error'
21 call error(3223)    !   stop 'place error'
22 call error(3224)    !   stop 'sqf error'
23 call error(3225)    !   stop 'restart error'
24 call error(3226)    !   stop 'vtk error'
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

subroutine construct(fnames,conc,place,snum,anglestempl,atempl,vallist,bead_type,ntlength,sphrad,graftdens)
! #################################################################
! #                                                               #
! #    subroutine 34:                                             #
! #    construct initial system                                   #
! #                                                               #
! #################################################################  
implicit none

common /cella/ dlxa,dlya,dlza
common /den/ rho

real*4,pointer :: rxt(:),ryt(:),rzt(:)
integer*1,pointer :: kindpt(:)
integer*4,pointer :: b1(:),b2(:),bt(:)
integer*4,pointer :: cn(:,:),cn2(:,:)

real*4 rho
real*4 conc(*)
character*16 fnames(*)
integer*1 snum,atempl,place(*)

integer*1 anglestempl(1000,3),vallist(10)
character*2 bead_type(*)

character*3 kindptc
character*256 str
character*52 dummy
integer trimtext
integer*4 npallt,natmsall,nbondsall,nat,nbd,nmol,nbcur,nbmax,nbmax2,nanglesall
integer*4 i,j,k,l,i1,i2,i3
real*4 dlxa,dlya,dlza
real*4 uni,rndposx,rndposy,rndposz
real*4 graftdens
integer*4 ntlength
real*4 sphrad,rsq,rb
logical exists

! use 3 temporary files; then combine them to the final file
open(1, file = 'atoms.dat',status='replace')
open(2, file = 'bonds.dat',status='replace')
open(3, file = 'angles.dat',status='replace') 

! calculate maximum number of beads in the box
npallt=dlxa*dlya*dlza*rho
natmsall=0
nbondsall=0
nanglesall=0
nbcur=0

! loop over species files
do k=1,snum
    if (fnames(k)(1:trimtext(fnames(k)))=='nanotube.gen') then

        nat=ntlength*7-1
        nbd=25*ntlength-20
        nmol=npallt*conc(k)/nat/100.0
    
        allocate (rxt(nat),ryt(nat),rzt(nat),kindpt(nat))
        allocate (b1(nbd),b2(nbd),bt(nbd))

        do j=0,ntlength-1
            rxt(j*6+1)=0
            ryt(j*6+1)=1.0
            rzt(j*6+1)=real(j)

            rxt(j*6+2)=sqrt(3.0)/2
            ryt(j*6+2)=0.5
            rzt(j*6+2)=real(j)

            rxt(j*6+3)=sqrt(3.0)/2
            ryt(j*6+3)=-0.5
            rzt(j*6+3)=real(j)

            rxt(j*6+4)=0
            ryt(j*6+4)=-1.0
            rzt(j*6+4)=real(j)

            rxt(j*6+5)=-sqrt(3.0)/2
            ryt(j*6+5)=-0.5
            rzt(j*6+5)=real(j)

            rxt(j*6+6)=-sqrt(3.0)/2
            ryt(j*6+6)=0.5
            rzt(j*6+6)=real(j)
    
        end do
        do j=0,ntlength-2
            rxt(ntlength*6+1+j)=0
            ryt(ntlength*6+1+j)=0
            rzt(ntlength*6+1+j)=real(j)+0.5
        end do
        kindpt=6
        kindpt(1:6*ntlength)=3

        !bonds
        do j=0,ntlength-1
            do i=1, 5
                b1(j*6+i)=i+j*6
                b2(j*6+i)=i+1+j*6
            end do
            b1(j*6+6)=6+j*6
            b2(j*6+6)=1+j*6
        end do


        do j=0,ntlength-2
            do i=1, 6
                b1(ntlength*6+i+j*6)=i+j*6
                b2(ntlength*6+i+j*6)=i+6+j*6
            end do

        end do


        do i=1, ntlength-2
            b1(ntlength*12-6+i)=ntlength*6+i
            b2(ntlength*12-6+i)=ntlength*6+1+i
        end do

        do j=0,ntlength-2
            do i=1,12
                b1(ntlength*13-8+i+j*12)=i+j*6
                b2(ntlength*13-8+i+j*12)=ntlength*6+1+j
            end do
        end do

        do i=1,nmol
        
            call rndrot(nat,rxt,ryt,rzt)
            rndposx=(uni())*dlxa
            rndposy=(uni())*dlya
            rndposz=(uni())*dlza
        
            do j=1,nat
                if (kindpt(j)==3) then
                    if (uni()<graftdens) then
                        write(1,'(2i4,3f14.6)')vallist(1), 1,rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                    else
                        write(1,'(2i4,3f14.6)')0, kindpt(j),rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                    end if
                else
                    write(1,'(2i4,3f14.6)')0, kindpt(j),rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                end if
                natmsall=natmsall+1
            end do

            ! output bonds
            do j=1,nbd
                write(2,'(2i8)'),b1(j)+nbcur,b2(j)+nbcur
                nbondsall=nbondsall+1
            end do
        
            do j=1,ntlength-2
                write(3,'(3i8)') (j-1)*6+1+nbcur,(j-1)*6+1+6+nbcur,(j-1)*6+1+12+nbcur
                write(3,'(3i8)') (j-1)*6+1+1+nbcur,(j-1)*6+1+7+nbcur,(j-1)*6+1+13+nbcur
                write(3,'(3i8)') (j-1)*6+1+2+nbcur,(j-1)*6+1+8+nbcur,(j-1)*6+1+14+nbcur
                write(3,'(3i8)') (j-1)*6+1+3+nbcur,(j-1)*6+1+9+nbcur,(j-1)*6+1+15+nbcur
                write(3,'(3i8)') (j-1)*6+1+4+nbcur,(j-1)*6+1+10+nbcur,(j-1)*6+1+16+nbcur
                write(3,'(3i8)') (j-1)*6+1+5+nbcur,(j-1)*6+1+11+nbcur,(j-1)*6+1+17+nbcur
                nanglesall=nanglesall+6
            end do
            do j=1,ntlength-3
                write(3,'(3i8)') ntlength*6+j+nbcur,ntlength*6+j+1+nbcur,ntlength*6+j+2+nbcur
                nanglesall=nanglesall+1
            end do

            nbcur=nat+nbcur
        end do

        ! delete temporary arrays
        nullify(rxt,ryt,rzt,kindpt,b1,b2,bt)
    elseif (fnames(k)(1:trimtext(fnames(k)))=='sphere.gen') then
        rb=0.7
        nat=0
        do i=-int(sphrad/0.61)-1,int(sphrad/0.61)+1
            do j=-int(sphrad/0.61)-1,int(sphrad/0.61)+1
                do l=-int(sphrad/0.61)-1,int(sphrad/0.61)+1
                    rsq=sqrt((i*0.61)**2+(j*0.61)**2+(l*0.61)**2)
                    if (rsq>sphrad) cycle    
                    nat=nat+1
                end do
            end do
        end do

        nmol=npallt*conc(k)/nat/100.0

        allocate (rxt(nat),ryt(nat),rzt(nat),kindpt(nat))
        kindpt=3
        nat=0

        do i=-int(sphrad/0.61)-1,int(sphrad/0.61)+1
            do j=-int(sphrad/0.61)-1,int(sphrad/0.61)+1
                do l=-int(sphrad/0.61)-1,int(sphrad/0.61)+1
                    rsq=sqrt((i*0.61)**2+(j*0.61)**2+(l*0.61)**2)
                    if (rsq>sphrad) cycle    
                    nat=nat+1
                    rxt(nat)=i*0.61
                    ryt(nat)=j*0.61
                    rzt(nat)=l*0.61
                end do
            end do
        end do

    
        nbd=0
        do i=1,nat
            do j=i+1,nat
                rsq=sqrt((rxt(i)-rxt(j))**2+(ryt(i)-ryt(j))**2+(rzt(i)-rzt(j))**2)
                if (rsq<=rb) nbd=nbd+1
            end do
        end do
    
        allocate (b1(nbd),b2(nbd),bt(nbd))
        allocate(cn(nat,0:0))
        cn=0
        nbd=0
        do i=1,nat
            do j=i+1,nat
                rsq=sqrt((rxt(i)-rxt(j))**2+(ryt(i)-ryt(j))**2+(rzt(i)-rzt(j))**2)
                if (rsq<=rb) then
                    nbd=nbd+1
                    b1(nbd)=i
                    b2(nbd)=j
                    cn(i,0)=cn(i,0)+1
                    cn(j,0)=cn(j,0)+1
                end if
            end do
        end do

        do i=1,nmol
            call rndrot(nat,rxt,ryt,rzt)
            rndposx=(uni())*dlxa
            rndposy=(uni())*dlya
            rndposz=(uni())*dlza
        
            do j=1,nat
        
                ! функционализация по кол-ву связей
                
                if (cn(j,0)<6) then
                    if (uni()<graftdens) then
                        write(1,'(2i4,3f14.6)')vallist(1), 1,rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                    else
                        write(1,'(2i4,3f14.6)')0, kindpt(j),rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                    end if
                end if
                if (cn(j,0)==6) write(1,'(2i4,3f14.6)')0,6,rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                natmsall=natmsall+1
            end do

            ! output bonds
            do j=1,nbd
                write(2,'(2i8)'),b1(j)+nbcur,b2(j)+nbcur
                nbondsall=nbondsall+1
            end do

            nbcur=nat+nbcur
        end do
        nullify(rxt,ryt,rzt,kindpt,b1,b2,bt,cn)

    else
        inquire( file = fnames(k)(1:trimtext(fnames(k))), exist = exists )
        if ( .not. exists ) call error(341)
        open(10, file = fnames(k)(1:trimtext(fnames(k))))

        ! read basic information about current species
        do
            read( 10, '(a)') str
            if ( str(40:44)=='V2000' ) then
                read(str(1:8),*) nat
                read(str(9:16),*) nbd
                exit
            end if
        end do
    
        ! calculate number of molecules of current species in the cell
        nmol=npallt*conc(k)/100/nat
    
        ! allocate arrays
        allocate (rxt(nat),ryt(nat),rzt(nat),kindpt(nat))
        allocate (b1(nbd),b2(nbd),bt(nbd))

        ! read bead structure
        do i = 1, nat
            read(10,'(3f10.4,a3)') rxt(i),ryt(i),rzt(i),kindptc
            do j=1,10
                if (kindptc(2:3)==bead_type(j)) exit
            end do
            if (j==11) call error(340)
            kindpt(i)=j
        end do 
        allocate(cn(nat,0:0),cn2(nat,0:0)) 
        cn=0
        cn2=0

    ! read bond structure
        do i = 1, nbd
            read(10,*)b1(i),b2(i),bt(i)
            if(bt(i)==2)cn(b1(i),0)=cn(b1(i),0)+1
            if(bt(i)==2)cn(b2(i),0)=cn(b2(i),0)+1
            cn2(b1(i),0)=cn2(b1(i),0)+1
            cn2(b2(i),0)=cn2(b2(i),0)+1
        end do
        nbmax=maxval(cn(:,0))
        nbmax2=maxval(cn2(:,0))
        nullify (cn,cn2)
        allocate (cn(nat,0:nbmax),cn2(nat,0:nbmax2))

        ! constuct connection matrixes for 2 types of bonds: for all bonds and only for double bonds
        cn=0
        cn2=0
        do i = 1, nbd
            i1=b1(i)
            i2=b2(i)
            if(bt(i)==2) then
              cn(i1,0) = cn(i1,0) + 1
              cn(i1,cn(i1,0)) = i2
              cn(i2,0) = cn(i2,0) + 1
              cn(i2,cn(i2,0)) = i1
            end if
            cn2(i1,0) = cn2(i1,0) + 1
            cn2(i1,cn2(i1,0)) = i2
            cn2(i2,0) = cn2(i2,0) + 1
            cn2(i2,cn2(i2,0)) = i1
        end do
    
        ! place molecules in the simulation cell
        do i=1,nmol
        
            ! random rotation variant 
            if (place(k)==1) then
                call rndrot(nat,rxt,ryt,rzt)
                rndposx=(uni())*dlxa
                rndposy=(uni())*dlya
                rndposz=(uni())*dlza

	    else if (place(k)==3) then
		!use coordinates from .mol files
		else if (place(k)==4) then
				call rndplace( nat,rxt,ryt,rzt,cn2,dlxa,dlya,dlza, 4 )
                rndposx = 0
                rndposy = 0
                rndposz = 0
            ! average distance between beads 0.5
            else
                call rndplace( nat,rxt,ryt,rzt,cn2,dlxa,dlya,dlza, 2 )
                rndposx = 0
                rndposy = 0
                rndposz = 0
            end if
        
            do j=1,nat
            
                ! output information about bead
                write(1,'(2i4,3f14.6)')vallist(kindpt(j)), kindpt(j),rxt(j)+rndposx,ryt(j)+rndposy,rzt(j)+rndposz
                natmsall=natmsall+1

                ! find angles types connected with this bead; output angles
                do i1=1,atempl
                    if (kindpt(j)==anglestempl(i1,1)) then
                        do i2=1,cn(j,0)
                            if (kindpt(cn(j,i2))==anglestempl(i1,2)) then
                                do i3=1, cn(cn(j,i2),0)
                                    if (kindpt(cn(cn(j,i2),i3))==anglestempl(i1,3)) then
                                        if (cn(cn(j,i2),i3)>j) then
                                            write(3,'(3i8)') j+nbcur,cn(j,i2)+nbcur,cn(cn(j,i2),i3)+nbcur
                                            nanglesall=nanglesall+1
                                        end if
                                    end if
                                end do
                            end if
                        end do
                        cycle
                    end if
                    if (kindpt(j)==anglestempl(i1,3)) then
                        do i2=1,cn(j,0)
                            if (kindpt(cn(j,i2))==anglestempl(i1,2)) then
                                do i3=1, cn(cn(j,i2),0)
                                    if (kindpt(cn(cn(j,i2),i3))==anglestempl(i1,1)) then
                                        if (cn(cn(j,i2),i3)>j) then
                                            write(3,'(3i8)') j+nbcur,cn(j,i2)+nbcur,cn(cn(j,i2),i3)+nbcur
                                            nanglesall=nanglesall+1
                                        end if
                                    end if
                                end do
                            end if
                        end do
                    end if
                end do
            end do

            ! output bonds
            do j=1,nbd
                write(2,'(2i8)'),b1(j)+nbcur,b2(j)+nbcur
                nbondsall=nbondsall+1
            end do
            nbcur=nat+nbcur
        end do

        ! delete temporary arrays
        nullify (cn,cn2)
        nullify(rxt,ryt,rzt,kindpt,b1,b2,bt)
    end if
end do
close(3)
close(2)
close(1)

! move information from temporary files to final file
open(1, file = 'atoms.dat')
open(2, file = 'bonds.dat') 
open(3, file = 'restart.dat',status='replace' )
open(4, file = 'angles.dat') 

! move general information
write(3,'(i8,f8.4)') natmsall, rho
write(3,'(3f20.12)') dlxa,dlya,dlza

! move information about beads
do i = 1, natmsall
    read (1,'(a50)') dummy
    write(3,'(i8,a50)')i,dummy
end do
close(1,status='delete')

! move information about bonds
write(3,'(a9,i9)') ' bonds:  ',nbondsall
do i=1,nbondsall
    read (2,'(a16)') dummy
    write(3,'(a16)') dummy
end do 
close(2,status='delete')

! move information about angles
write(3,'(a9,i9)') ' angles:  ',nanglesall
do i=1,nanglesall
    read (4,'(a24)') dummy
    write(3,'(a24)') dummy
end do 
close(4,status='delete')
close(3)
end

subroutine rndrot(nat,x,y,z)
! #################################################################
! #                                                               #
! #    subroutine 35:                                             #
! #    rotate molecule through random angle                       #
! #    around the centre of mass                                  #
! #                                                               #
! #################################################################  
dimension x(*),y(*),z(*)
dimension rot(9),qtn(4)
real*4 uni

!construct random quaternions

qtn(1)=(2.0*uni()-1.0)
qtn(2)=(2.0*uni()-1.0)
qtn(3)=(2.0*uni()-1.0)
qtn(4)=(2.0*uni()-1.0)

rnm=1.0/sqrt(qtn(1)**2+qtn(2)**2+qtn(3)**2+qtn(4)**2)

qtn(1)=rnm*qtn(1)
qtn(2)=rnm*qtn(2)
qtn(3)=rnm*qtn(3)
qtn(4)=rnm*qtn(4)

!construct rotation matrix

rot(1) = qtn(1)**2+qtn(2)**2-qtn(3)**2-qtn(4)**2
rot(2) = 2.d0*(qtn(2)*qtn(3) - qtn(1)*qtn(4))
rot(3) = 2.d0*(qtn(2)*qtn(4) + qtn(1)*qtn(3))
rot(4) = 2.d0*(qtn(2)*qtn(3) + qtn(1)*qtn(4))
rot(5) = qtn(1)**2-qtn(2)**2+qtn(3)**2-qtn(4)**2
rot(6) = 2.d0*(qtn(3)*qtn(4) - qtn(1)*qtn(2))
rot(7) = 2.d0*(qtn(2)*qtn(4) - qtn(1)*qtn(3))
rot(8) = 2.d0*(qtn(3)*qtn(4) + qtn(1)*qtn(2))
rot(9) = qtn(1)**2-qtn(2)**2-qtn(3)**2+qtn(4)**2

!calculate centre of mass

cma=0.0
cmx=0.0
cmy=0.0
cmz=0.0

do i=1,nat
    cma=cma+1.0
    cmx=cmx+x(i)
    cmy=cmy+y(i)
    cmz=cmz+z(i)
end do

cmx=cmx /cma
cmy=cmy /cma
cmz=cmz /cma

!rotate cluster about centre of mass

do i=1,nat
    xs=x(i)-cmx
    ys=y(i)-cmy
    zs=z(i)-cmz
    x(i)=rot(1)*xs+rot(4)*ys+rot(7)*zs+cmx
    y(i)=rot(2)*xs+rot(5)*ys+rot(8)*zs+cmy
    z(i)=rot(3)*xs+rot(6)*ys+rot(9)*zs+cmz
end do


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
do i = 1, nat

    ! exit if the bead has no bonds
    if ( cn(i,0) == 0 ) then
        id(i)=-1
        cycle
    end if

    ! place all connected to the ith bead beads next to it
    do k = 1, cn(i,0)
        j = cn(i,k)

        ! don't move jth bead if it has already been moved 
        if (id(j)==-1) cycle
        dx = ( uni() - 0.5 )
        dy = ( uni() - 0.5 )
        dz = ( uni() - 0.5 )
        r = sqrt( dx**2 + dy**2 + dz**2 )
        if ( r==0.0 ) then
		r = 1.0
		end if
		if (place==4) then
        scale = 0.5 / r
        dx = dx * scale
        dy = dy * scale
        dz = dz * scale
        rxt(j)=rxt(i)+dx
        ryt(j)=ryt(i)+dy
        rzt(j)=rzt(i)+dz
		else
		scale = 1.0 / r
        dx = dx * scale
        dy = dy * scale
        dz = dz * scale
        rxt(j)=rxt(i)+dx
        ryt(j)=ryt(i)+dy
        rzt(j)=rzt(i)+dz
		end if
		!write(*,*) dx, dy, dz, scale, r, place
        id(j) = -1
    end do
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
if( flag==120) then
    call barrier
    if(rank==0) stop 'error 120:too many VTK frames to output'
    stop
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
if( flag==190) then
    stop 'error 190:communication array overflow in forcesangle subroutine'
end if
if( flag==251) then
    stop 'error 251:domain size along y became less than Rcut during deformation'
end if
if( flag==252) then
    stop 'error 252:domain size along z became less than Rcut during deformation'
end if
if( flag==253) then
    stop 'error 253:domain size along x became less than Rcut during deformation'
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
if( flag==325) then
    stop 'error 325: input script error, sum of concentrations /= 100'
end if
if( flag==324) then
    stop 'error 324: input script error, number of concentrations /= number of species'
end if
if( flag==323) then
    stop 'error 323: input script error, number of species is equal to 0'
end if
if( flag==326) then
    stop 'error 326: input script error, species_number command error'
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
if( flag==3212) then
    stop 'error 3212: input script error, steps_3 command error'
end if
if( flag==3213) then
    stop 'error 3213: input script error, steps_ch command error'
end if
if( flag==3214) then
    stop 'error 3214: input script error, pair_coeff command error'
end if
if( flag==3215) then
    stop 'error 3215: input script error, chem_coeff command error'
end if
if( flag==3216) then
    stop 'error 3216: input script error, there is no dpdconf.dat file'
end if
if( flag==3217) then
    stop 'error 3217: input script error, density command error'
end if
if( flag==3218) then
    stop 'error 3218: input script error, output_freq command error'
end if
if( flag==3219) then
    stop 'error 3219: input script error, angle_template command error'
end if
if( flag==3220) then
    stop 'error 3219: input script error, angle_template command error'
end if
if( flag==3220) then
    stop 'error 3220: input script error, bond_template command error'
end if
if( flag==3221) then
    stop 'error 3221: input script error, valency command error'
end if
if( flag==3222) then
    stop 'error 3222: input script error, deform command error'
end if
if( flag==3223) then
    stop 'error 3223: input script error, place command error'
end if
if( flag==3224) then
    stop 'error 3224: input script error, sqf command error'
end if
if( flag==3225) then
    stop 'error 3225: input script error, restart command error'
end if
if( flag==3226) then
    stop 'error 3226: input script error, vtk command error'
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
if( flag==341) then
    stop 'error 341:input .mol file does not exist'
end if
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