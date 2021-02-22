using Plots





function show_analytic_model(f::Function)
    plot(-10:0.1:10, f.(-10:0.1:10), linewidth=5)
end

function analytic_dos(f::Function; mesh=100, histogram_width=100, max_energy=10)

    dos_array = zeros(histogram_width*max_energy)
    for x_iter in 1:mesh
        for y_iter in 1:mesh
            x = x_iter/mesh*100
            y = y_iter/mesh*100
            energy = f(x, y)

            if energy<max_energy
                dos_array[round(Int, energy*histogram_width)] = dos_array[round(Int, energy*histogram_width)] + histogram_width/mesh^2
            end

        end
    end

    return dos_array

end