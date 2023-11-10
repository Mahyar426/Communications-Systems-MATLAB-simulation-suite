function [dist_profile,dist,Ai]=DistanceProfile(k,crc_encoder)
v=zeros(2^k,k);
for i=1:(2^k)-1
    v(i+1,:)=de2bi(i,k);
    c(i+1,:)=(crc_encoder(v(i+1,:)'))';
end
for i=1:size(c,1)
    wh(i,:)=c(1,:)+c(i,:);
end
wh=mod(wh,2);
dist_profile=sum(wh,2);
dist=min(dist_profile((2:end),:));
Ai=zeros(1+2^k,1);

for i=1:size(dist_profile,1)
    Ai(dist_profile(i,1)+1,1)=Ai(dist_profile(i,1)+1,1)+1;
end

end