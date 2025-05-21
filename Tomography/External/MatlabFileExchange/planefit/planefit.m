function varargout=planefit(xyz)
%Fit 3D plane to given (x,y,z) data
%
%    [a,b,c,d] =  planefit(xyz)   
%     [abc,d]  =  planefit(xyz)
%      hcoeff  =  planefit(xyz)
%           
%IN:        
%           
%    xyz:   3xN input array with coordinates organized in columns
%           [x1,x2,x3...;y1,y2,y3,...;z1,z2,z3,...]
%                 
%OUT:       
%           
%    [a,b,c,d] : Coefficients of fitted plane equation a*x+b*y+c*z=d
%    [abc,d]   : Coefficients in 2-argument form where abc=[a,b,c]
%     hcoeff   : Homogeneous coefficient vector hcoeff=[a,b,c,-d]


 if size(xyz,1)~=3
    error 'Input xyz matrix must be 3xN' 
 end

 xyz=xyz.';
 
 mu=mean(xyz,1);
 
 [~,~,V]=svd(xyz-mu,0);
 
 normal=V(:,end).';
 d=normal*mu';
 
 switch nargout
     
     case {0,1}
 
      varargout={[normal,-d]};
      
     case 2
         
      varargout={ normal, d};
      
     case 4
       
      varargout=[num2cell(normal),{d}]; 
         
     otherwise
         
      error 'Must be 1,2, or 4 output args'
          
 end
 
end