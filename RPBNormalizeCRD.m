function [crd] = RPBNormalizeCRD(crd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  if max(crd(1,:)) ~= 0
        crd(1,:)=0.5*crd(1,:)./(max(crd(1,:)));
    end
    if max(crd(2,:)) ~= 0
        crd(2,:)=0.5*crd(2,:)./(max(crd(2,:)));
    end
    if max(crd(3,:)) ~= 0
        crd(3,:)=0.5*crd(3,:)./(max(crd(3,:)));
    end
    if min(crd(1,:))~= 0
        crd(1,:)=-0.5*crd(1,:)./(min(crd(1,:)));
    end
    if min(crd(2,:)) ~= 0
        crd(2,:)=-0.5*crd(2,:)./(min(crd(2,:)));
    end
    if min(crd(3,:)) ~= 0
       crd(3,:)=-0.5*crd(3,:)./(min(crd(3,:)));
    end
end

