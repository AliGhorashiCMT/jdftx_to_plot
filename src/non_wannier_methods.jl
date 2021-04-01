function nonwannier3dimepsilon(filebase::String, lattice::Array{<:Array{<:Real, 1}, 1},  numkpoints::Integer, numbands::Integer; spin::Integer=1, histogram_width::Real=10)
    
    momentumfile = "$filebase.momenta"
    eigfile = "$filebase.eigenvals"
    momenta = np.reshape(np.fromfile(momentumfile, dtype=np.complex), (numkpoints, 3, numbands, numbands))
    momentasquared = sum(abs.(momenta), dims=2)[:, 1, :, : ] ##Get rid of extra axis now that summation has been performed
    energies = np.reshape(np.fromfile(eigfile), (numkpoints, numbands)) ./ eV
    prefactor = e²ϵ*π*(spin/3)*ħ^2/(mₑ^2) #Factor of three takes into account isotropy. Spin is conventionally taken to be 2
    Epsilons = zeros(100*histogram_width)
    Vol = unit_cell_volume(lattice_vectors)

    for k in 1:numkpoints
        for band1 in 1:numbands
            for band2 in 1:numbands
                energy1, energy2 = energies[k, band1], energies[k, band2]
                ω  = energy2-energy1
                if ω>0
                    pabs = momentasquared[k, band1, band2]*(ħ/bohrtoangstrom)^2
                    f2 = heaviside(μ-energy2)
                    f1 = heaviside(μ-energy1)
                    Epsilons[round(Int, histogram_width*ω)+1] = Epsilons[round(Int, histogram_width*ω)+1] + (f1-f2)*prefactor*(1/Vol)*pabs*histogram_width*(1/ω^2)*1/numkpoints
                end
            end
        end
    end
    return Epsilons
end

function nonwannierimpol(filebase::String, lattice::Array{<:Array{<:Real, 1}, 1},  q::Real, numkpoints::Integer, numbands::Integer, ::Val{3}; spin::Integer=1, histogram_width::Real=10)
    momentumfile = "$filebase.momenta"
    eigfile = "$filebase.eigenvals"

    momenta = np.reshape(np.fromfile(momentumfile, dtype=np.complex), (numkpoints, 3, numbands, numbands))
    momentasquared = sum(abs.(momenta), dims=2)[:, 1, :, : ] ##Get rid of extra axis now that summation has been performed
    
    energies = np.reshape(np.fromfile(eigfile), (numkpoints, numbands)) ./ eV
    
    Impols = zeros(100*histogram_width)
    V = unit_cell_volume(lattice_vectors)

    for k in 1:numkpoints
        for band1 in 1:numbands
            for band2 in 1:numbands

                band1 == band2 && continue ##Don't consider intraband transitions

                energy1, energy2 = energies[k, band1], energies[k, band2]
                ω  = energy2-energy1
                if ω>0
                    pabs = momentasquared[k, band1, band2]
                    f2 = heaviside(μ-energy2)
                    f1 = heaviside(μ-energy1)

                    overlap = q^2ħ^4/bohrtoangstrom^2*1/(mₑ)^2*1/(ω^2)

                    Impols[round(Int, histogram_width*ω)+1] = Impols[round(Int, histogram_width*ω)+1] + π*(f2-f1)/V*overlap*(1/numkpoints)*histogram_width*spin

                end
            end
        end
    end
    return Impols
end


function nonwannierimpol(filebase::String, lattice::Array{<:Array{<:Real, 1}, 1},  q::Real, numkpoints::Integer, numbands::Integer, ::Val{2}; spin::Integer=1, histogram_width::Real=10)
    momentumfile = "$filebase.momenta"
    eigfile = "$filebase.eigenvals"

    momenta = np.reshape(np.fromfile(momentumfile, dtype=np.complex), (numkpoints, 3, numbands, numbands))
    momentasquared = sum(abs.(momenta), dims=2)[:, 1, :, : ] ##Get rid of extra axis now that summation has been performed
    
    energies = np.reshape(np.fromfile(eigfile), (numkpoints, numbands)) ./ eV
    
    Impols = zeros(100*histogram_width)
    V = unit_cell_area(lattice_vectors)

    for k in 1:numkpoints
        for band1 in 1:numbands
            for band2 in 1:numbands

                band1 == band2 && continue ##Don't consider intraband transitions

                energy1, energy2 = energies[k, band1], energies[k, band2]
                ω  = energy2-energy1
                if ω>0
                    pabs = momentasquared[k, band1, band2]
                    f2 = heaviside(μ-energy2)
                    f1 = heaviside(μ-energy1)

                    overlap = q^2ħ^4/bohrtoangstrom^2*1/(mₑ)^2*1/(ω^2)

                    Impols[round(Int, histogram_width*ω)+1] = Impols[round(Int, histogram_width*ω)+1] + π*(f2-f1)/V*overlap*(1/numkpoints)*histogram_width*spin

                end
            end
        end
    end
    return Impols
ends