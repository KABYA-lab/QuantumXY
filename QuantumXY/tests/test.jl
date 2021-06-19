module QuantumXY





using Test
export QuantumXY
using LinearAlgebra





@testset "GS_h_<0" begin
    HH=Hamiltonian_basis(N=3,J=0,h=-1.0)
    evalues,evectors=eigen(HH)
    @test evectors[:,1] ==  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] 
end


@testset "GS_h_>0" begin
    HH=Hamiltonian_basis(N=3,J=0,h=1.0)
    evalues,evectors=eigen(HH)
    @test evectors[:,1] == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end


end
