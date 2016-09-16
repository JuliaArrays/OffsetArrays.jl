subroutine setprob

    implicit none
    character*25 :: fname
    integer :: iunit
    real(kind=8) :: rho,bulk,cc,zz,beta

    common /cparam/ rho,bulk,cc,zz   
    common /cqinit/ beta
 
    ! Set the material parameters for the acoustic equations
    ! Passed to the Riemann solver rp1.f in a common block
 
    iunit = 7
    fname = 'setprob.data'
    ! open the unit with new routine from Clawpack 4.4 to skip over
    ! comment lines starting with #:
    call opendatafile(iunit, fname)


    ! density:
    read(7,*) rho

    ! bulk modulus:
    read(7,*) bulk

    ! sound speed:
    cc = dsqrt(bulk/rho)

    ! impedance:
    zz = cc*rho

    ! beta for initial conditions:
    read(7,*) beta

end subroutine setprob
