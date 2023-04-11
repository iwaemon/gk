program main
    implicit none
    integer::ier,Nt,i,j,sn
    double precision::eps
    complex(kind(0d0))::a,b,s,s0
    complex(kind(0d0))::xa(1:5),xb(1:5)
    double precision :: t 
    
    eps=1d-10
    ier=1
  
    Nt=1000
    do i=1,Nt
       t=i*(1d0-(-1d0))/dble(Nt)-1d0
       s=dcmplx(0d0,0d0)
       
       sn=1
       if(t.gt.0d0)sn=-1
       if(t.eq.0)sn=0
       
       xa(1)=dcmplx(-2000d0,sn*2000d0)
       xb(1)=dcmplx(-1d0,0d0)
  
       xa(2)=dcmplx(-1d0,0d0)
       xb(2)=dcmplx(-1d0,1d0)
  
       xa(3)=dcmplx(-1d0,1d0)
       xb(3)=dcmplx( 1d0,1d0)
  
       xa(4)=dcmplx( 1d0,1d0)
       xb(4)=dcmplx( 1d0,0d0)
  
       xa(5)=dcmplx( 1d0,0d0)
       xb(5)=dcmplx( 2000d0,sn*2000d0)  
  
       do j=1,5
          s0=dcmplx(0d0,0d0)
          call cqag_sk(wrapper,xa(j),xb(j),eps,s0,3,ier)
          s=s+s0
       enddo
       write(10,*)t,dble(s),dimag(s)
    enddo  
       
    stop
    contains
    function wrapper(z)
        implicit none
        complex(kind(0d0)):: wrapper
        complex(kind(0d0)),intent(in)::z
        complex(kind(0d0)), external :: g 
        wrapper = g(z,t)
      end function wrapper
  end program main
  
  function g(z,t)
    implicit none
    complex(kind(0d0))::g
    complex(kind(0d0)),intent(in)::z
    double precision,intent(in)::t
  
    double precision::pi=dacos(-1d0)
   
    g=exp(-dcmplx(0d0,1d0)*z*t)/(z)
    g=g*(-1d0/(2d0*pi*dcmplx(0d0,1d0)))
  
    return
  end function g