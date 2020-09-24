module random
  
  implicit none

  ! Benny: trick for easy debugging (reproducing results with a constant seed)
  integer, save :: seed = 0

  contains

  subroutine set_seed(new_seed)
  
    implicit none
    
    integer, intent(IN) :: new_seed
    
    seed = new_seed
  
  end subroutine set_seed
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function randold()

    implicit none

    integer :: time

    real :: randold

    ! Benny: debugging
    !print *, "seed", seed    
    if (seed.eq.0) seed = time()

    randold = ran3(seed)

  end function randold

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Added by Benny, include a pid argument
  ! From now 05102013 replaces the above rand(old) with the following rand
  ! After implementing: it seems all processors are referring to the same seed.
  function rand()

    implicit none

    integer :: time, getpid

    real :: rand

    !print *, seed
!    if (seed.eq.0) seed = time()
    !if (seed.eq.0) seed = ABS( MOD((time() * 181) * ((getpid() - 83) * 359), 104729))
    !if (seed.eq.0) seed = ABS(getpid())
!    rand = ran3(seed)

    if (seed .eq. 0) then
      seed = time()
      call random_seed()
    end if

    call random_number(rand)

  end function rand

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function randnozero()

    implicit none

    integer :: time, getpid

    real :: randnozero
    
    if (seed .eq. 0) then
      seed = time()
      call random_seed()
    end if

    do 
      call random_number(randnozero)
      if (randnozero .ne. 0.0d0) return
    end do

  end function randnozero

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Numerical Recipes ran3, updated for Fortran 90

  FUNCTION ran3(idum)  
        INTEGER :: idum  
        INTEGER :: MBIG,MSEED,MZ  
  !C     REAL :: MBIG,MSEED,MZ  
        REAL :: ran3,FAC  
        PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)  
  !C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)  
        INTEGER :: i,iff,ii,inext,inextp,k  
        INTEGER :: mj,mk,ma(55)  
  !C     REAL :: mj,mk,ma(55)  
        SAVE iff,inext,inextp,ma  
        DATA iff /0/  
        if(idum.lt.0.or.iff.eq.0)then  
          iff=1  
          ! Benny: there seems to be potential bug reported for ran3 here
          ! iabs is the integer absolute value function.
          ! mj = iabs(MSEED-iabs(idum)) ! make sure mj is not negative
          mj=MSEED-iabs(idum)  
          mj=mod(mj,MBIG)  
          ma(55)=mj  
          mk=1  
          do 11 i=1,54  
            ii=mod(21*i,55)  
            ma(ii)=mk  
            mk=mj-mk  
            if(mk.lt.MZ)mk=mk+MBIG  
            mj=ma(ii)  
  11      continue  
          do 13 k=1,4  
            do 12 i=1,55  
              ma(i)=ma(i)-ma(1+mod(i+30,55))  
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG  
  12        continue  
  13      continue  
          inext=0  
          inextp=31  
          idum=1  
        endif  
        inext=inext+1  
        if(inext.eq.56)inext=1  
        inextp=inextp+1  
        if(inextp.eq.56)inextp=1  
        mj=ma(inext)-ma(inextp)  
        if(mj.lt.MZ)mj=mj+MBIG  
        ma(inext)=mj  
        ran3=mj*FAC  
        return  
  END FUNCTION ran3

end module random
