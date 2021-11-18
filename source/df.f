      function df(e,zr,ar,zl,alm break)                                  
c     ******************************************************************  
c     damage function using the lindhard partition of                    
c     energy between atomic and electronic motion                        
c     call with e=0 for each reaction to precompute the constants      
c     ******************************************************************                                                  
      implicit real*8 (a-h,o-z)  
      real e, zr, ar, zl,alm, break, df                                           
      data twothd/.666666667d0/                                           
      data threeq/.75d0/                                                  
      data sixth/.166666667d0/                                          
      data onep5/1.5d0/                                                  
      data c1/30.724d0/                                                   
      data c2/.0793d0/                                                    
      data c3/3.4008d0/                                                   
      data c4/.40244d0/                                                  
      zero=0                                                              
c                                                                        
c                                                                           
c     Normal displacement energy algorithm from Kinchin-Pease            
c                                                                         
      if (zr.eq.0.) go to 120                                         
      if (e.gt.0.) go to 100                                            
      el=c1*zr*zl*sqrt(zr**twothd+zl**twothd)*(ar+al)/al                
      rel=1/el                                                           
      fl=c2*zr**twothd*sqrt(zl)*(ar+al)**onep5/                        
     1   ((zr**twothd+zl**twothd)**threeq*ar**onep5*sqrt(al))            
      df=0.                                                              
      return                                                          
  100 if (e.lt.break) go to 120                                       
      ep=e*rel                                                       
      dam=e/(1+fl*(c3*ep**sixth+c4*ep**threeq+ep))                    
      df=dam                                                        
      return                                                            
  120 df=0.                                                             
      return                                                             
      end                                                                
