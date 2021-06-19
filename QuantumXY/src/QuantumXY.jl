module QuantumXY




export bit_repr, generate_basis, Hamiltonian_basis, Hamiltonian_kron, magnetization, correlation



using LinearAlgebra





function bit_repr(n::Integer, N::Integer)
    aa= BitArray(parse(Bool, i) for i in string(n, base=2, pad=N))
    return aa
end


       

 
#generates the basis vector        
function generate_basis(N::Integer)
    nstates = 2^N
    basis = Vector{BitArray{1}}(undef, nstates)
    for i in 0:nstates-1
        basis[i+1] = bit_repr(i, N)
    end
    return reverse(basis)
end






function Hamiltonian_basis(;N::Integer,J::Integer,h::Float64)
    basis=generate_basis(N)
    H0 = zeros(Float64,2^N,2^N)
    for (i,bstate) in enumerate(basis)
        for j in 1:N
            if j<N
                new_bstate=copy(bstate)
                if new_bstate[j]!=new_bstate[j+1]
                    new_bstate[j],new_bstate[j+1]=new_bstate[j+1],new_bstate[j]
                    new_i = findfirst(x->x == new_bstate,basis)
                    H0[i,new_i] -= J
                elseif new_bstate[j]==new_bstate[j+1]
                    new_i = findfirst(x->x == new_bstate,basis)
                    H0[i,new_i] =0.0
                end
            elseif j==N
                new_bstate=copy(bstate)
                if new_bstate[j]!=new_bstate[1]
                    new_bstate[j],new_bstate[1]=new_bstate[1],new_bstate[j]
                    new_i = findfirst(x->x == new_bstate,basis)
                    H0[i,new_i] -= J
                elseif new_bstate[j]==new_bstate[1]
                    new_i = findfirst(x->x == new_bstate,basis)
                    H0[i,new_i] =0.0
                end
            end
        end
    end
    
    HB = zeros(Float64,2^N,2^N)
    for i in 1:2^N
        HB[i,i] -= h*(count(x->x==1,basis[i])-count(x->x==0,basis[i]))
    end
    return (H0+HB)
end







function Hamiltonian_kron(;N::Integer,J::Integer,h::Float64)
    II = [1 0; 0 1]
    σp = [0 1; 0 0]
    σm = [0 0; 1 0]
    σz = [1 0; 0 -1]
    

    fst_t = fill(II, N)
    fst_t[1] = σp
    fst_t[2] = σm
    
    
    fst_tc = fill(II,N)
    fst_tc[1] = σm
    fst_tc[2] = σp
    

    snd_t = fill(II, N)
    snd_t[1] = σz
    
    H = zeros(Int, 2^N, 2^N)
    for i in 1:N
        H -= J*foldl(kron, fst_t)
        fst_t = circshift(fst_t,1)
    end
    

    for i in 1:N
        H -= J*foldl(kron, fst_tc)
        fst_tc = circshift(fst_tc,1)
    end
    
    
    for i in 1:N
        H -= h*foldl(kron, snd_t)
        snd_t = circshift(snd_t,1)
    end
    H
end






function magnetization(state, basis)
    M = 0.
    for (i, bstate) in enumerate(basis)
        bstate_M = 0.
        for spin in bstate
            bstate_M += (state[i]^2 * (spin ? 1 : -1))/length(bstate)
        end
        M += abs(bstate_M)
    end
    return M
end





function correlation(state, basis,r)
    CC = 0
    for (i, bstate) in enumerate(basis)
        bstate_CC = 0
        for (s,s1) in zip(bstate, circshift(bstate,r))
            bstate_CC += (state[i]^2 * (s ? 1 : -1)*(s1 ? 1 : -1))/length(bstate)
        end
        CC += bstate_CC
    end
    MM = magnetization(state,basis)    
    
    return CC-(MM*MM)
end


end
















