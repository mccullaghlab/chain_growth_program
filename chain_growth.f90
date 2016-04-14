
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Modules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ompData

        integer nThreads

endmodule ompData

module chainData

        integer nMonomers
        real (kind=8) bond_k
        real (kind=8) bond_eq
        real (kind=8) ang_k
        real (kind=8) ang_eq
        real (kind=8) dih_k
        real (kind=8) dih_eq
        real (kind=8) sphere_r 
        real (kind=8) sphere_k
        real (kind=8), allocatable :: chainCoord(:,:,:)

endmodule chainData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program chain_growth
        use chainData
        use ompData
        implicit none
        character*80 outFile, parmFile
        real (kind=8) omp_get_wtime
        real (kind=8) ti,tf

        ti = omp_get_wtime()

        call parse_command_line(nMonomers,sphere_r,parmFile,outFile,nThreads)

        call read_parm_file(parmFile)

        call grow_the_chain(outFile)

        tf = omp_get_wtime()
        write(*,'("Total time elapsed:",f12.4)') tf-ti

endprogram chain_growth



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine grow_the_chain(outFile)
        use chainData
        use ompData
        implicit none
        integer, parameter :: nChains = 100 
        integer, parameter :: nPoints = 60
        real (kind=8) newMonomerCoord(nPoints,3)
        real (kind=8) newMonomerEnergy(nPoints)
        character*80 outFile
        real (kind=8) r, theta, phi
        real (kind=8) spherical(3), temp(3)
        integer monomer, i, j, k, points
        integer minPoint
        real (kind=8) minEnergy
        real (kind=8) chainEnergy(nChains)
        integer chain, minChain

        allocate(chainCoord(nChains,nMonomers,3))
        chainEnergy = 0
        !$omp parallel private(r,chain, points,spherical,temp,monomer,newMonomerCoord,newMonomerEnergy,minEnergy,minPoint) shared(nMonomers,chainEnergy,bond_eq,chainCoord) num_threads(nThreads)
        !$omp do 
        do chain = 1, nChains
                write(*,'("Working on chain: ", i5)') chain
                flush(6)
               ! pick random starting location for first bead
               call random_theta_phi(spherical(3), spherical(2))
               call random_number(r)
               spherical(1) = sphere_r*r
               call spherical_to_cartesian(spherical, chainCoord(chain,1,:))
               ! pick random direction for second bead
               call random_theta_phi(spherical(3), spherical(2))
               spherical(1) = bond_eq
               call spherical_to_cartesian(spherical, temp)
               chainCoord(chain,2,:) = chainCoord(chain,1,:) + temp
               
               do monomer = 3, nMonomers
        
                       do points = 1, nPoints
        
                               call fibonacci_sphere(nPoints,points, temp, bond_eq)
                               newMonomerCoord(points,:) = chainCoord(chain,monomer-1,:) + temp
        
                               call compute_new_monomer_energy(chainCoord(chain,1:monomer,:),monomer,newMonomerCoord(points,:),newMonomerEnergy(points))
        
                       enddo
        
                       do points = 1, nPoints
                               if (points==1 .or. newMonomerEnergy(points) < minEnergy) then
                                       minEnergy = newMonomerEnergy(points)
                                       minPoint = points
                               endif
                       enddo
                       chainEnergy(chain) = chainEnergy(chain) + minEnergy
                       chainCoord(chain,monomer,:) = newMonomerCoord(minPoint,:)
                enddo 
        enddo
        !$omp end do nowait
        !$omp end parallel
        
        do chain=1,nChains
                print*, chain, chainEnergy(chain)
                if (chain==1 .or. chainEnergy(chain) < minEnergy) then
                        minEnergy = chainEnergy(chain)
                        minChain = chain
                endif
        enddo
        print*, "minChain = ", minChain
        open(24,file=outFile)
        write(24,'(i10)') nMonomers
        write(24,'(i10)') nMonomers
        do monomer=1,nMonomers
                write(24,'("C  ",4f12.4)') chainCoord(minChain,monomer,1), chainCoord(minChain,monomer,2),chainCoord(minChain,monomer,3),sqrt(dot_product(chainCoord(minChain,monomer,:),chainCoord(minChain,monomer,:)))
        enddo
        close(24)

endsubroutine grow_the_chain

subroutine fibonacci_sphere(nPoints, point, xyz, r_seed)
        implicit none
        real (kind=8), parameter :: pi=3.1415926535
        integer nPoints
        integer point
        real (kind=8) xyz(3)
        real (kind=8) offset
        real (kind=8) increment
        real (kind=8) r
        real (kind=8) r_seed
        real (kind=8) phi

        
        offset = 2.0/dble(nPoints)
        increment = pi * (3.0-sqrt(5.0))

        xyz(2) = (point-1)*offset - 1 + offset/2.0

        r = sqrt(1-xyz(2)**2)

        phi = (point-1) * increment

        xyz(1) = cos(phi) * r
        xyz(3) = sin(phi) * r

        xyz = xyz * r_seed

endsubroutine fibonacci_sphere

subroutine compute_new_monomer_energy(chainCoord,nMonomers,newCoord,newEnergy)
        use chainData, only : sphere_r
        implicit none
!        real (kind=8), parameter :: a = 5e-15
!        real (kind=8), parameter :: a = 5e-14
        real (kind=8), parameter :: a = 10.0
!        real (kind=8), parameter :: sigma6 = 244140625.0
!        real (kind=8), parameter :: sigma12 = 5.9604645e16
!        real (kind=8), parameter :: sigma6 = 1771561.0
!        real (kind=8), parameter :: sigma12 = 3138428376721.0
        real (kind=8), parameter :: charge_factor = 16.6025  ! dielectric of 80.0
!        real (kind=8), parameter :: charge_factor = 25.000
!        real (kind=8), parameter :: sigma6 = 262144.0
!        real (kind=8), parameter :: sigma12 = 68719476736.0
        real (kind=8), parameter :: eps = 1.00
        real (kind=8), parameter :: sigma = 25.0
        real (kind=8) sigma12
        integer nMonomers
        real (kind=8) chainCoord(nMonomers,3)
        real (kind=8) newCoord(3)
        real (kind=8) newEnergy
        real (kind=8) r2, r6, r12
        real (kind=8) diff(3)
        integer monomer

        sigma12 = sigma**12

        newEnergy = 0

        !angle energy
        call compute_angle_energy(newCoord,chainCoord(nMonomers,:),chainCoord(nMonomers-1,:),newEnergy)

        ! pair energy
        if (nMonomers>4) then
                do monomer = 1, nMonomers-4 
                        diff = chainCoord(monomer,:) - newCoord
                        r2 = dot_product(diff,diff)
                        if (monomer < nMonomers-10) then
                                r6 = r2*r2*r2
                                r12 = r6*r6
                                newEnergy = newEnergy + eps *  sigma12/r12 
                        endif
                        newEnergy = newEnergy +  charge_factor/sqrt(r2)

                enddo
        endif

        ! now add spherical boundary potential
        r2 = dot_product(newCoord,newCoord)
        ! to the sixth power
!        newEnergy = newEnergy + a * r2*r2*r2
        newEnergy = newEnergy + (a / (sphere_r-sqrt(r2))) ** 12

endsubroutine

! compute angle energy
subroutine compute_angle_energy(r1,r2,r3,energy)
        use chainData
        implicit none
        real (kind=8) r1(3), r2(3), r3(3)
        real (kind=8) r21(3), r23(3)
        real (kind=8) energy
        real (kind=8) theta

        r21 = r1 - r2
        r21 = r21 / sqrt(dot_product(r21,r21))
        r23 = r3 - r2
        r23 = r23 / sqrt(dot_product(r23,r23))
        theta = acos(dot_product(r21,r23))

        energy = ang_k * (theta - ang_eq)**2
        
endsubroutine compute_angle_energy

! generate uniformly distributed theta and phi values
subroutine random_theta_phi(theta,phi)
        implicit none
        real (kind=8), parameter :: pi = 3.1415926535
        real (kind=8) u, v,theta, phi

        call random_number(u)
        theta = 2 * pi * u
        call random_number(v)
        phi = acos(2 * v - 1)

endsubroutine random_theta_phi

!convert from spherical to cartesian
subroutine spherical_to_cartesian(spherical, cartesian)
        implicit none
        real (kind=8) spherical(3), cartesian(3)

        cartesian(1) = spherical(1) * sin(spherical(2)) * cos(spherical(3))
        cartesian(2) = spherical(1) * sin(spherical(2)) * sin(spherical(3))
        cartesian(3) = spherical(1) * cos(spherical(2)) 

endsubroutine spherical_to_cartesian

!read command line information
subroutine parse_command_line(nMonomers,sphere_r,parmFile,outFile,nThreads)
        implicit none
        integer i
        character*30 arg
        character*80 parmFile
        character*80 outFile
        integer nMonomers
        integer nThreads
        real (kind=8) sphere_r
        logical nMonomerFlag
        logical sphereRFlag
        logical parmFileFlag
        logical outFileFlag
        logical nThreadsFlag

        nMonomerFlag = .false.
        sphereRFlag = .false.
        outFileFlag = .false.
        parmFileFlag = .false.
        nThreadsFlag = .false.
        i=1
        do 
   
                call get_command_argument(i, arg) 
   
                select case (arg) 
   
                case ('-parm')
                        i = i+1
                        call get_command_argument(i,parmFile)
                        parmFileFlag=.true.
                        print*, "PARM File: ", parmFile
                case ('-o')
                        i = i+1
                        call get_command_argument(i,outFile)
                        outFileFlag=.true.
                        print*, "Output data file: ", outFile
                case ('-nm') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') nMonomers
                        print*, "Number of monomers to be grown:", nMonomers
                        nMonomerFlag = .true.
                case ('-r') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(f10.5)') sphere_r
                        print*, "Radius of sphere:", sphere_r
                        sphereRFlag = .true.
                case ('-np') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') nThreads
                        print*, "Number of threads: ", nThreads
                        nThreadsFlag = .true.
                case default 
                        print '(a,a,/)', 'Unrecognized command-line option: ', arg 
                        print*, 'Usage: chain_growth.x -parm [parm file] -nm [number of monomers] -o [output file] -r [sphere radius] -np [number of threads]'
                        stop 
                end select 
                i = i+1
                if (i.ge.command_argument_count()) exit
        enddo

        if (parmFileFlag.eqv..false.) then
                write(*,'("No parameter file was given. Must provide parameter file using command line argument -parm [parm file name]")')
                stop
        endif
        if (nMonomerFlag.eqv..false.) then
                write(*,'("Number of monomers was not give.  Must provide number of monomers using command line argument -nm [number of monomers]")')
                stop
        endif
        if (outFileFlag.eqv..false.) then
                write(*,'("Must provide an output file name  using command line argument -o [output file name]")')
                stop
        endif
        if (sphereRFlag.eqv..false.) then
                write(*,'("Radius of sphere to grow chain in was not given.  Must provide using command line argument -r [radius of sphere]")')
                stop
        endif
        if (nThreadsFlag.eqv..false.) then
                write(*,'("Using default of 1 thread.  Change this with command line argument -np [number of threads]")')
                nThreads = 1
        endif

endsubroutine parse_command_line


!Read atomic charges from some file
subroutine read_parm_file(parmFile)
        use chainData
        implicit none
        integer atom, j
        character*6 check           !character to check if NATOM is in the line
        character*8 numChar         !character to read number of atoms.  must be converted to an integer
        character*4 atomCheck
        character*24 posChar
        character*80 parmFile
        integer ios

        !open the psf file
        open(10,file=parmFile)


        close(10)

        ang_k = 0.00332
        ang_eq = 2.76
        dih_k = 0.00332
        dih_eq = 2.76
        bond_eq = 2.8

endsubroutine read_parm_file




