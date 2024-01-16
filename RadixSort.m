function Asorted=RadixSort(A,digit)
i=digit;
Asorted=A;
while i>=1
    [x,inx]=sort(Asorted(:,i));
    Asorted(:,:)=Asorted(inx(:),:);
    i=i-1;
end
end %function