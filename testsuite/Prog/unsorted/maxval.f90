Program test1

Complex(Kind = 8) :: arr(100)
Integer :: i
Real(Kind = 8) :: tmp

Do I=1,100
arr(i) = CMPLX(i,i,kind=8)
enddo
tmp = 1000.0
DO I = 1,100
           if (abs(dble(arr(I))) <   tmp ) tmp = abs(dble(arr(I)))
        ENDDO
write (*,*) tmp

tmp = minval (abs(dble(arr)))
write (*,*) tmp

end program test1
