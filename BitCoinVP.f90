!****************************************************************************
!                                                                           *
!  PROGRAM: BitCoinValuePredictor                                           *
!                                                                           *
!  PURPOSE:  This program, will try to predict the future value of a        *
!  cryptovalue, by using a statistical approach.                            *
!  In particular, using a .dat file containing more than 2100 hystorical    *
!  value for BitCoin, it will define some usefull statistical variables,    *
!  by taking n value where n is choosed from the user.                      *
!  The defined variables, will be used to predict the future value of n+1   *
!  bitcoin value, and will be compared to the real one.                     *
!  The predicted values and the real values will then be compared using a   *
!  Ki^2 test, to define the best n for the predictions, and see how well    *
!  the approach will esteem the future value.                               *
!                                                                           *
!****************************************************************************

! Main Program

    program BitCoinVP

    implicit none

    ! Variables !

    integer (kind=4) :: n,i,correct
    integer (kind=4), parameter :: NMAX = 2900

    character (len=1) :: c

    real (kind=8) :: xavg,sigx,nsigx,xpred,kisq,TSTART,TSTOP,correctperc

    real(kind=8),dimension(:),allocatable :: x,xapp

    logical (kind=1) opt

    ! Body of BitCoinVP !

    ! Defining Initial Value !

    opt= .TRUE.

    do while(opt)

        i = 0
        correct = 0
        xavg = 0.
        sigx = 0.
        kisq = 0.
        n = NMAX + 100

        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

        ! Defining the number of bitcoin values used for the statistical analysis

        do while(n <= 0 .or. n > NMAX)
           
            write (*,*) 'Please insert the number of values to be used in the analysis (between 0 and ',NMAX,' :'
    	    read (*,*) n

        end do 
        
        ! Defining the dimension of dynamical vectors

        allocate (x(n),xapp(n))

        ! Opening the Bitcoin datafile

        open (unit = 1, file = 'bitcoindata.txt')

        ! Defining the reading format from txt file

101     format(f9.2)

        ! Opening the result datafile

        open (unit = 2, file = 'BitSimulationResults.txt')
    	write (2,*) '################################################################'
   	    write (2,*) '#              Results of the numerical simulation           #'
    	write (2,*) '################################################################'
    	write (2,*) ''
    	write (2,*) '#   record n°          truevalue            simulatedvalue     #'
    	write (2,*) ''

    	!Defining the writing format on results file 

201   format (i5,a,f9.2,a,f9.2)
        
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

        ! Initial CPU Time

        call CPU_TIME(TSTART)

        write (*,*) 'The simulation has started'

        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
 
        ! Defining the values vector and the actual average value

        do while (i < n)
        
            i = i + 1
        	read(1,101) x(i)
        	xavg = xavg + x(i)

        end do

        xavg = xavg/n
       
        ! Defining the remaining statistical variables

        call InvertArray(n,x,xapp)

        call SigmaCalc(n,xavg,x,sigx)

        nsigx = (x(1) - xavg)/(sigx**(1./2.))

        xpred =  x(1) + ((nsigx*xavg)/10.)

       ! Beginning of the prediction cycle

       i = 0

        do while (i < (NMAX - n))
       
              i = i + 1
              read(1,101) xapp(1)
              write (2,201) i,' ',xapp(1),' ',xpred
              kisq = kisq + (((xapp(1) - xpred)**2)/xpred)
              xavg = xavg - (x(n)/n) + (xapp(1)/n)

              if ((nsigx < 0 .and. (xapp(1) - x(1)) <=0 ) .or. (nsigx > 0 .and. (xapp(1) - x(1)) >=0 )) then

              	  correct = correct + 1

              endif

              call ShiftArray (n,x,xapp)

              call SigmaCalc(n,xavg,x,sigx)

              nsigx = (x(1) - xavg)/(sigx**(1./2.))

              xpred =  x(1) + ((nsigx*xavg)/20.)


        enddo

        kisq = kisq/(NMAX - n)

       correctperc = ((correct*100.)/(NMAX-n))

       ! Deallocating the dynamical vectors

       deallocate (x,xapp)

       ! Closing files
       close (unit = 1)
       close (unit = 2)
       
       ! Ending CPU Time

       call CPU_TIME (TSTOP)

       ! Execution Completed!

        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Simulation completed in ',TSTOP-TSTART,' seconds'
        write (*,*) 'The value obtained of Ki squared is : ',kisq
        write (*,*) 'The software predicted ',correct,' values on ',NMAX-n,' total values, resulting in a '&
        ,correctperc,' percent correct predictions'

        ! Relaunching request

999   	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Do you wish to relaunch the program (y/n)?'
        read (*,*) c
    	if(c/='y'.and.c/='n'.and.c/='Y'.and.c/='N') then
		    write (*,*) 'Wrong Immission!'
            goto 999
   	    endif
        
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        if(c=='n'.or.c=='N') then
		        opt = .FALSE.
   	    endif
    
    end do

    ! End Program !

    end program

    ! SUBROUTINE !

    !This subroutine, invert a given n dimensional array using a auxiliary one 

    subroutine InvertArray (n,x,xapp)

        integer (kind=4) :: n,i
        real(kind=8),dimension(n) :: x,xapp

        i = 0
        do while (i < (n/2))
         
             i = i + 1
             xapp(n - i + 1) = x(i)
             x(i) = x(n - i + 1)
             x(n - i + 1) = xapp(n - i + 1)

        enddo

        return

    end

    ! This subroutine, given an array of values and the average value, return the value of sigma squared

    subroutine SigmaCalc(n,xavg,x,sigma)

        integer (kind=4) :: n,i
        real(kind=8) :: xavg,sigma
        real(kind=8),dimension(n) :: x

        i = 0
        sigma = 0.

        do while (i < n)
        
            i = i + 1
            sigma = sigma + (x(i) - xavg)**2  

        end do

        sigma = sigma/n

        return

    end

    ! This subroutine, memorize the first value of the auxiliary array on the first slot of the main array, then shift all the other positions of the main array by one

    subroutine ShiftArray (n,x,xapp)

        integer (kind=4) :: n,i
        real(kind=8),dimension(n) :: x,xapp

        i = 0
        do while (i < (n-1))
         
             i = i + 1
             xapp(i + 1) = x(i)
             x(i) = xapp(i)

        enddo

        x(n) = xapp(n) 

        return

    end