function Listnew=generatesisd(Listold)

%This function implements the procedure of Bracho et al to generate all
%three dimensional strongly involutive self dual (sisd) polytopes with a given
%number of vertices given the list of all such polytopes with one less vertex. 

%    [Javier Bracho, Luis Montejano, Eric Prez, and Jorge Ramrez Alfonsn.
%     Strongly involutive self-dual polyhedra. 
%     Ars Mathematica Contemporanea, 20(1):143{149, 2021.]

% Observation: The most unefficient part of this algorithm is the procedure
% to avoid redundancy in the generated polytopes. This currently relies on a graph
% isomorphism check, using the generic matlab's isisomorphic function. This should be
% easy to improve with a more specific procedure tailored to the current
% problem, but this version was more than enough to generate a sufficient
% number of examples.


%Input: 
%
%Listold - a cell array whose entries are the symmetric supports of all sisd
%3-polytopes with a fixed number of vertices
%
%Output: 
%
%Listnew - a cell array whose entries are the symmetric supports of all sisd
%3-polytopes with a number of vertices one more than those of Listold

nmat=length(Listold);
n=length(Listold{1});

Listnew={};
Splitv={};
t=0;

%For every matrix in Listold, generate all possible splits along facets
%with more than 4 vertices

for index=1:nmat
    A=Listold{index};
    for i=1:n   %Running over all facets
        
        I=find(~A(i,:)); %I is the unordered set of vertices in the facet
        
        k=length(I);
        if k>=4  %If the facet has at least 4 vertices
            
            %Find J the set of facets that intersect with facet i along edges 
            J=zeros(k);
            auxi=0;
            for j=1:n
                if sum(A(j,I)==0)==2
                    auxi=auxi+1;
                    J(auxi)=j;
                end
            end
            
            %Now we have to reorder I and J into I2 and J2 so that adjacent
            %indices correspond to adjacent vertices/facets
            
            I2=[I(1)];
            J2=[];
            for j=1:k
                if A(I(1),J(j))==0
                    J2(1)=J(j);
                    break;
                end
            end
            for j=2:k
                for l=1:k
                    if A(I(l),J2(j-1))==0&&(~(I(l)==I2(j-1)))                        
                        I2(j)=I(l);
                        break;
                    end
                end
                
                for l=1:k
                    if A(I2(j),J(l))==0&&(~(J(l)==J2(j-1)))
                        J2(j)=J(l);
                        break;
                    end
                end
            end
            
            %Now that we have this, we are ready to split the facet i along
            %every possible cord of the cycle.
            
            
            A2=[A(1:i,:);A(i:n,:)];
            A2=[A2(:,1:i), A2(:,i:n)];	
            for j=1:k-2   
                for l=j+2:k-2+min([j,2]) %Iterate over all cords (j,l) of the k-cycle
                    A3=A2;
                    %We start by generating a new candidate slack matrix A3 by splitting the cycle with new edge at (j,l)
                    for ind=1:j-1
                        aux=I2(ind);
                        if aux>i
                            aux=aux+1;
                        end
                        A3(i+1,aux)=1;
                        A3(aux,i+1)=1;
                    end
                    for ind=j+1:l-1
                        aux=I2(ind);
                        if aux>i
                            aux=aux+1;
                        end
                        A3(i,aux)=1;
                        A3(aux,i)=1;
                    end
                    for ind=l+1:k
                        aux=I2(ind);
                        if aux>i
                            aux=aux+1;
                        end
                        A3(i+1,aux)=1;
                        A3(aux,i+1)=1;
                    end
                    
                    %We now have to check if the new candidate is already
                    %on the list or not
                    
                    %v will be the ordered vector of the sizes of the
                    %facets of the candidate new polytope. This will allow
                    %a preliminary check, to avoid too many isomorphism
                    %tests.
                    
                    v=zeros(n+1,1);
                    for ind=1:n+1
                        for ind2=1:n+1
                            if A3(ind,ind2)==0
                                v(ind)=v(ind)+1;
                            end
                        end
                    end
                    v=sort(v);   
                    
                    %We genertae the graph associated to the candidate
                    %matrix and check if it is isomorphic to any of the
                    %graphs of the previously generated matrices
                    G=graph([zeros(n+1) A3;A3' zeros(n+1)]);
                    addmat=true;
                    for ind=1:t
                        if isequal(v,Splitv{ind})
                            if isisomorphic(G,graph([zeros(n+1) Listnew{ind};Listnew{ind}' zeros(n+1)]))
                                addmat=false;
                                break;
                            end
                        end
                    end
                    %If it is not, then we add it to the list
                    if addmat
                        t=t+1;
                        Listnew{t}=A3;
                        Splitv{t}=v;
                    end
                end
            end
        end
    end
end

%If the polytopes to be listed have an even number of vertices then we have
%to add the pyramid over a polygon to the list, as it cannot be obtained by
%splitting
if mod(n+1,2)==0
    t=t+1;
    Listnew{t}=[1, zeros(1,n); zeros(n,1), ones(n)-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1)-diag(1,n-1)-diag(1,-n+1)];
end

