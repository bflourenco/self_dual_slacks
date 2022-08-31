function [issisd,Msym]=issisd(M)

%This function receives the support of a slack matrix of a (combinatorial) polytope. It
%checks if said polytope is strongly involutively self dual (sisd), and if so
%gives back a symmetric matrix Msym, with ones on the diagonal, attained by
%permuting the columns of M and that certifies that M is sisd.


n=size(M,1);
Msym=zeros(n);


%STEP 1 - As an heuristical speed up we start by listing all feasible possibilities
%for the last two columns. We need them to have the same numbers of ones
%has the last two rows and ones on the two last diagonal positions, as well
%as making the entries (n,n-1) and (n-1,n) equal.
%The possible pairs of indices for these two last columns will be stored as
%rows of a matrix vv.
aux=sum(M(n,:));
v=[];
for j=1:n
      if M(n,j)==1&&sum(M(:,j))==aux
          v=[v,j];
      end
end
ncand=length(v);

vv=[];
aux=sum(M(n-1,:));
for i=1:ncand
    for j=1:n
        if M(n-1,j)==1&&sum(M(:,j))==aux&&M(n-1,v(i))==M(n,j)&&~(v(i)==j)
           vv=[vv;j v(i)]; 
        end
    end
end
n2cand=size(vv,1);


%STEP 2 - for each feasible pair found in STEP 1 we will generate all
%candidates for indices of each of the remaining columns. 
%
%Column j must have the same number of ones as the row j, and have one on position j.
%They also must make the submatrix indexed by {j,n-1,n} symmetric.
%
%For each j we generate that list of candidates, and store it has
%candidatos{j,i}, where i is the index of the choice for the last two
%entries that we are using. w2 is a vector that stores the size of all
%these lists.

candidatos=cell(n,n2cand);
w2=zeros(n-2,n2cand);

for i=1:n2cand
    for j=1:n-2
        aux=sum(M(j,:));
        v=[];
        for k=1:n
            if M(j,k)==1&&sum(M(:,k))==aux&&M(j,vv(i,2))==M(n,k)&&M(j,vv(i,1))==M(n-1,k)&&~(k==vv(i,2))&&~(k==vv(i,1))
                v=[v,k];
            end
        end
        if isempty(v)
            break;
        else
            candidatos{j,i}=v;
            w2(j,i)=length(v);
        end
    end
    candidatos{n-1,i}=[vv(i,1)];
    candidatos{n,i}=[vv(i,2)];
end

%Step 3 - Now we remove any pair that resulted in some vertex having no
%possible candidates.

count=sum(min(w2)>0);
k=1;
candidatos3=cell(n,count);
for i=1:count
   while min(w2(:,k))==0
       k=k+1;
   end
   for j=1:n
       candidatos3(j,i)=candidatos(j,k);
   end
   k=k+1;
end


for ii=1:count
    %For every choice of the last two entries, we have, for each row index, 
    %a list of vectors of possible column indices that could be put in that
    %position. We now want to generate all possible simultaneous choices of
    %distinct indices for these columns. I.e., we want a permutation p such
    %that entry p_i is in the candidate list for vertex i.
    
    %In order to do this we will implement a vectorial counter that will go 
    %through all possible choices and check if they verify our property.
    
    auxi=zeros(n,1);
    w=zeros(n-1,1);
    for i=1:n
        auxi(i)=candidatos3{i,ii}(1);
        w(i)=length(candidatos3{i,ii});
    end
    
    %Inicialize counter
    v = ones(1, n);  % Index vector
    
    impossible=0;
    
    while impossible==0
        
        k=n-1;
        % Search for repetitions from right to left
        while k>0 %while there are repetitions
            if ismember(auxi(k),auxi(k+1:n)) %if k is repeated
             v(k)=v(k)+1; %advance candidate
             if v(k) > w(k) %there are no more candidates
                 %Increase counter
                 v(k)=1;
                 for j=k+1:n
                     v(j)=v(j)+1;
                     if v(j)<=w(j)
                         auxi(j)=candidatos3{j,ii}(v(j));
                         k=j;
                         break;
                     end
                     if v(j)>w(j) && j~=n
                         v(j)=1;
                         k=j;
                         auxi(j)=candidatos3{j,ii}(v(j));
                     end
                 end
             end
             
             %If we advanced k, we reset all previous entries to 1
             for j=1:k-1
                 v(j)=1;
                 auxi(j)=candidatos3{j,ii}(1);
             end
             auxi(k)=candidatos3{k,ii}(v(k));
             
            else
                k=k-1;
            end
            %If the counter is spent we exit
            if v(n)>w(n)
                impossible=1;
                break;
            end
        end
        
        if impossible==1
            break;
        end
        
        %There was a feasible permutation
        %Then we check if the corresponding matrix is symmetric, and if it
        %is we stop our search
        
        Maux=M(:,auxi(:));
        if max(max(abs(Maux-Maux')))==0 
            Msym=Maux;
            issisd=true;
            break;
        end
        
        %If it is not symmetric we increase the counter and resume the
        %search
        v(1)=v(1)+1;
        if v(1) > w(1)
            v(1)=1;
            for l=2:n
                v(l)=v(l)+1;
                if v(l)<=w(l)
                    auxi(l)=candidatos3{l,ii}(v(l));
                    break;
                end
                if v(l)>w(l) && l~=n
                    v(l)=1;
                    auxi(l)=candidatos3{l,ii}(v(l));
                end
            end
        end
        auxi(1)=candidatos3{1,ii}(v(1));
    end
    %If we already have what we want we break the search
    if issisd
        break;
    end
end








