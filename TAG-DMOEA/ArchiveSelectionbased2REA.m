function [Archive1,Archive2,W,EI,znad,B] = ArchiveSelectionbased2REA(Population1,Population2,W,N,EI,znad,B,WA,RA,z)
  [FrontNO,~]       = NDSort(Population1.objs,N);
  Archive1          = Population1(FrontNO==1);
  if RA==1
      [DistanceValue,~,PopObj1]  = F_distance(Archive1.objs);                                  
      [~,rank]                = sort(DistanceValue,'ascend');    
      if length(Archive1)>N
       Archive1               = PopObj1(rank(1:N),:);
      else
       Archive1               = PopObj1; 
      end
      znew              = max(Archive1,[],1);
      znad              = (0.3*znad+0.7*znew);
       Archive1         = [];
  end
  if WA==1 && length(find(EI>10))>ceil(N/10) 
      [FrontNO,~]        = NDSort(Population2.objs,N);
       Archive2          = Population2(FrontNO==1);
      if length(Archive2)>10*N 
          p                  = Shape_Estimate(Archive2,N,z,znad);
          [Archive2]         = EnvironmentalSelection2REA(Archive2,N,p,z,znad);
          W                  = Archive2.objs;
          EI                 = zeros(1,N);
          z                  = min(W,[],1);
          znad               = max(W,[],1);
          W                  = (W-repmat(z,N,1))./(repmat(znad-z,N,1));
          Archive2           = [];
            T                = ceil(N/10);
            B                = pdist2(W,W);
            [~,B]            = sort(B,2);
            B                = B(:,1:T);
      end
  else
        Archive2=[];
  end    
end

    