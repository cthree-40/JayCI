program test
        use detci2
        implicit none
        integer :: Address, Num_Elec, Num_Orb, Num_Det
        integer, dimension(:), allocatable :: Array

        Address = 52
        Num_Elec= 2
        Num_Orb = 15
        Num_Det = 105

        allocate( Array(2) )

        call genorbstring( Address, Num_Elec, Num_Orb, Num_Det, Array )

        print *, Array

end program test
