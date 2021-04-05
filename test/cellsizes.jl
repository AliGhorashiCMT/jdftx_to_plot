@testset "cell volumes" begin
    ##A trivial cubic lattice
    @test unit_cell_volume([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ 1
    ##An FCC lattice like silicon
    @test unit_cell_volume([[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]]) ≈ 1/4
    ##A BCC lattice
    @test unit_cell_volume([[1/2, 1/2, 1/2], [1/2, -1/2, 1/2], [1/2, 1/2, -1/2]]) ≈ 1/2
end

@testset "cell areas" begin
    ##A trivial square lattice
    @test unit_cell_area([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ 1
    a = 1.42*sqrt(3) ## graphene lattice parameter
    @test unit_cell_area([[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]) ≈ a^2*sin(π/3)  
    ##Test that giving distances in Bohr also works. Note that lattice automatically makes everything into Bohr for ease of JDFTX outputs
    ##And unit_cell_area always gives everything in Angstroms^2, so we convert back from Bohr to Anstroms. 
    @test(unit_cell_area(lattice([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))) ≈ 1
    @test unit_cell_area(lattice([[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]])) ≈ a^2*sin(π/3)  
    @test(unit_cell_area(lattice([1/bohrtoangstrom 0 0; 0 1/bohrtoangstrom 0; 0 0 1/bohrtoangstrom]))) ≈ 1
    @test unit_cell_area(lattice([a/bohrtoangstrom -a/2/bohrtoangstrom 0; 0 a*sqrt(3)/2/bohrtoangstrom  0; 0 0 10/bohrtoangstrom])) ≈ a^2*sin(π/3) 
end

@testset "cell brillouin zone areas" begin
    @test brillouin_zone_area([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) ≈ 4*π^2
    ##Brillouin zone area of graphene
    a = 1.42*sqrt(3) ## graphene lattice parameter
    K = 4π/(3a)
    @test brillouin_zone_area([[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]) ≈ 3*sqrt(3)/2*K^2
end