
!     Power 4
!     integer, parameter     :: L1_SIZE   = 32*1024
!     integer, parameter     :: L1_WAY    = 2
!     integer, parameter     :: L2_SIZE   = 512*1024
!     integer, parameter     :: L2_WAY    = 8
      
!     PENTIUM 4(Northwood or prior)
!     integer, parameter     :: L1_SIZE   = 8*1024
!     integer, parameter     :: L1_WAY    = 4
!     integer, parameter     :: L2_SIZE   = 512*1024
!     integer, parameter     :: L2_WAY    = 8
      
!     PENTIUM 4(Prescott)
!     integer, parameter     :: L1_SIZE   = 16*1024
!     integer, parameter     :: L1_WAY    = 8
!     integer, parameter     :: L2_SIZE   = 1024*1024
!     integer, parameter     :: L2_WAY    = 8
      
!     CERELON D(Prescott)
!     integer, parameter     :: L1_SIZE   = 16*1024
!     integer, parameter     :: L1_WAY    = 8
!     integer, parameter     :: L2_SIZE   = 256*1024
!     integer, parameter     :: L2_WAY    = 4

!     Itanium 2
!     integer, parameter     :: L1_SIZE   = 16*1024
!     integer, parameter     :: L1_WAY    = 8
!     integer, parameter     :: L2_SIZE   = 256*1024
!     integer, parameter     :: L2_WAY    = 8

!     Core2
!      integer, parameter     :: L1_SIZE   = 32*1024
!      integer, parameter     :: L1_WAY    = 8
!      integer, parameter     :: L2_SIZE   = 2*1024*1024
!      integer, parameter     :: L2_WAY    = 16

!     A64FX
      integer, parameter     :: L1_SIZE      = 64*1024
      integer, parameter     :: L1_WAY       = 4
      integer, parameter     :: L1_LINE_SIZE = 256
      integer, parameter     :: L2_SIZE      = 8*1024*1024
      integer, parameter     :: L2_WAY       = 16
      integer, parameter     :: L2_LINE_SIZE = 256


      integer, parameter     :: L1_LSIZE     = (L1_SIZE/L1_WAY)/8
      integer, parameter     :: L1_WINDOW    = L1_LINE_SIZE/8
      integer, parameter     :: L2_LSIZE     = (L2_SIZE/L2_WAY)/8
      integer, parameter     :: L2_WINDOW    = 2*L2_LINE_SIZE/8

      integer, parameter     :: n_columns = L2_LSIZE 
      integer, parameter     :: PREFETCH_SIZE = 512/8

!     fr IA32-Linux
      integer, parameter     :: PAGE_SIZE  = 4096
      integer, parameter     :: PAGE_LSIZE = PAGE_SIZE/8

