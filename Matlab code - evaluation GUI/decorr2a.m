function [du,dv,dc]=decorr2a(dm0,dm1)
[M,N]=size(dm0);
dm0=double(dm0); dm1=double(dm1);
c=abs(fftshift(ifft2(conj(fft2(dm0)).*fft2(dm1))));
cn=abs(fftshift(ifft2(conj(fft2(dm0))./abs(fft2(dm0)).*fft2(dm1)./abs(fft2(dm1)))));
   cc=sum(sum((dm0-dm1).*(dm0-dm1)));
   if cc==0 | isnan(sum(c(:))),
   	du=nan;
   	dv=nan;
    dc=nan;
   else 
      [dv,du]=find(c==max(c(:)));

      if dv==1|dv==M|du==1|du==N,
            	du=nan;
                dv=nan;
                dc=nan;
      else
        
        cc=c(dv-1:dv+1,du-1:du+1);
        [mm,nn]=size(cc);
        [xx,yy]=meshgrid(-(mm-1)/2:(mm-1)/2,-(nn-1)/2:(nn-1)/2);
        xx2=xx.*xx;
        yy2=yy.*yy;
        A= [sum(sum(xx2.*xx2)) sum(sum(xx2.*yy2)) 0 0 sum(sum(xx2));
            sum(sum(xx2.*yy2)) sum(sum(yy2.*yy2)) 0 0 sum(sum(yy2));
            0 0 sum(sum(xx2)) 0 0;
            0 0 0 sum(sum(yy2)) 0;
            sum(sum(xx2)) sum(sum(yy2)) 0 0 sum(sum(ones(mm,nn)))];
        y=[sum(sum(xx2.*cc)) sum(sum(yy2.*cc)) sum(sum(xx.*cc)) sum(sum(yy.*cc)) sum(sum(cc))]';
        x=A\y;
        x0=-x(3)/2/x(1);
        y0=-x(4)/2/x(2);
        dc=sqrt(x(1)^2+x(2)^2);

      du=du+x0-(N/2+1);
      dv=dv+y0-(M/2+1);
    end 
  end
