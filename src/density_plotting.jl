function plot_density(density_file::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}})
    n = np.fromfile(density_file, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
            break ##Only look at first instance of fftbox 
        end
    end
    V=0
    for r in readlines(outfile)
        if contains(r, "volume")
            splittedV = split(r)
            for split in splittedV
                try 
                    V = parse(Float64, split)
                catch

                end
            end
        end
    end
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    println("Nelectrons = ", dV * np.sum(n))
    n = np.reshape(n, S)
    if perpaxis isa Val{'x'}
        n = n[1,:,:]        
        n = np.roll(n, Int(S[2]/2), axis=0) 
        n = np.roll(n, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        n = n[:,1,:]        
        n = np.roll(n, Int(S[1]/2), axis=0) 
        n = np.roll(n, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        n = n[:,:,1]        
        n = np.roll(n, Int(S[1]/2), axis=0) 
        n = np.roll(n, Int(S[2]/2), axis=1) 
    end
    heatmap(n)
end

function plot_diffdensity(density_file1::String, density_file2::String, outfile::String, perpaxis::Union{Val{'x'}, Val{'y'}, Val{'z'}})
    nup = np.fromfile(density_file1, dtype=np.float64)   
    ndn = np.fromfile(density_file2, dtype=np.float64)   
    ##Obtain volume and Chosen FFT Box from output file
    S = Vector{Int}()
    for r in readlines(outfile)
        if contains(r, "Chosen fftbox")
            splittedfft = split(r)
            for split in splittedfft
                try
                    a = parse(Int, split)
                    push!(S, a)
                catch

                end
            end
            break ##Only look at first instance of fftbox 
        end
    end
    V=0
    for r in readlines(outfile)
        if contains(r, "volume")
            splittedV = split(r)
            for split in splittedV
                try 
                    V = parse(Float64, split)
                catch

                end
            end
        end
    end
    println("Unit Cell Volume: ", V)
    dV = V / np.prod(S)
    println("Nelectrons of Species 1= ", dV * np.sum(nup))
    println("Nelectrons of Species 2= ", dV * np.sum(ndn))
    nup = np.reshape(nup, S)
    ndn = np.reshape(ndn, S)
    if perpaxis isa Val{'x'}
        nup = nup[1,:,:]        
        nup = np.roll(nup, Int(S[2]/2), axis=0) 
        nup = np.roll(nup, Int(S[3]/2), axis=1) 
        ndn = ndn[1,:,:]        
        ndn = np.roll(ndn, Int(S[2]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'y'}
        nup = nup[:,1,:]        
        nup = np.roll(nup, Int(S[1]/2), axis=0) 
        nup = np.roll(nup, Int(S[3]/2), axis=1)
        ndn = ndn[:,1,:]        
        ndn = np.roll(ndn, Int(S[1]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[3]/2), axis=1) 
    elseif perpaxis isa Val{'z'}
        nup = nup[:,:,1]        
        nup = np.roll(nup, Int(S[1]/2), axis=0) 
        nup = np.roll(nup, Int(S[2]/2), axis=1) 
        ndn = ndn[:,:,1]        
        ndn = np.roll(ndn, Int(S[1]/2), axis=0) 
        ndn = np.roll(ndn, Int(S[2]/2), axis=1) 
    end
    heatmap(nup-ndn)
end