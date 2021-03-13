# We will include models for the dielectric function of graphene for future reference 

#We start first with the density of states of graphene


#=

The functions below are exact expressions for the dynamical polarization of graphene at charge neutrality 
and at finite doping. They are copied verbatim from: B Wunsch et al 2006 New J. Phys. 8 318


=#


function real_neutral(q, w)
    return -1*q^2/(4*(-w^2+36*q^2)^.5)
end

function imag_neutral(q, w)
    return -q^2/(4*(w^2-36*q^2)^.5)
end

function intraband_1a_real(q, w, mu)
    return -2*mu/(36*pi)+1/4*(q^2)/sqrt(abs(36*(q^2)-w^2))
end

function gplus(x)
    return x*sqrt(x^2-1)-log(x+sqrt(x^2-1))
end

function gminus(x)
    if x>=0
        return x*sqrt(1-x^2)-atan(sqrt(1-x^2)/x)
    else
        return -pi+x*sqrt(1-x^2)+atan(sqrt(1-x^2)/(-x))
    end
end

function f(q, w)
    return 1/(4*pi)*q^2/sqrt(abs(w^2-(6*q)^2))
end


function intraband_1b_real(q, w, mu)
    return -2*mu/(6^2*pi)+f(q, w)*(gplus((2*mu+w)/(6*q))-gplus((2*mu-w)/(6*q)))
end


function intraband_1b_imag(q, w, mu)
    return f(q, w)*pi
end


function intraband_1a_imag(q, w, mu)
    return f(q, w)*(gplus((2*mu-w)/(6*q))-gplus((2*mu+w)/(6*q)))
end


function intraband_2a_imag(q, w, mu)
    return -f(q, w)*gplus((2*mu+w)/(6*q))
end


function intraband_2b_imag(q, w, mu)
    return -f(q, w)*gminus((w-2*mu)/(6*q))
end

function intraband_2a_real(q, w, mu)
    return -2*mu/(36*pi)-f(q, w)*gminus((w-2*mu)/(6*q))
end

function intraband_2b_real(q, w, mu)
    return -2*mu/(36*pi)+f(q, w)*gplus((w+2*mu)/(6*q))
end

function intraband_3a_real(q, w, mu)
    return -2*mu/(36*pi)+f(q, w)*(gminus((2*mu+w)/(6*q))-gminus((w-2*mu)/(6*q)))
end

function intraband_3b_real(q, w, mu)
    return -2*mu/(36*pi)+f(q, w)*(gplus((2*mu+w)/(6*q))-gplus((w-2*mu)/(6*q)))
end

function intraband_real_total(q, w, mu)
    if w<6*q && w<2*mu-6*q
        return intraband_1a_real(q, w, mu)
    elseif w<-2*mu+6*q
        return intraband_3a_real(q, w, mu)
    elseif w<6*q && w>2*mu-6*q && w>-2*mu+6*q
        return intraband_2a_real(q, w, mu)
    elseif w>6*q && w<2*mu-6*q
        return intraband_1b_real(q, w, mu)
    elseif w>6*q && w>2*mu-6*q && w<2*mu+6*q
        return intraband_2b_real(q, w, mu)
    elseif w>2*mu+6*q
        return intraband_3b_real(q, w, mu)
    else
        return 0
    end
end

function intraband_imag_total(q, w, mu)
    if w<6*q && w<2*mu-6*q
        return intraband_1a_real(q, w, mu)
    elseif w<-2*mu+6*q
        return 0
    elseif w<6*q && w>2*mu-6*q && w>-2*mu+6*q
        return intraband_2a_imag(q, w, mu)
    elseif w>6*q && w<2*mu-6*q
        return intraband_1b_imag(q, w, mu)
    elseif w>6*q && w>2*mu-6*q && w<2*mu+6*q
        return intraband_2b_imag(q, w, mu)
    elseif w>2*mu+6*q
        return 0
    else
        return 0
    end
end

function graphene_total_polarization(q, w, mu)

    return intraband_real_total(q, w, mu)+ (w<6q ? real_neutral(q, w) : 0 )

end

function exact_graphene_epsilon(q, w, mu)
    return 1-e²ϵ/2/q*graphene_total_polarization(q, w, mu) 
end

function graphene_energy(t, kx, ky)
    t*sqrt(3+2*cos(sqrt(3)*kx*1.42)+4*cos(3/2*ky*1.42)*cos(sqrt(3)*kx/2*1.42))
end

function graphene_energy_normalizedk(t, graphene_lattice, k1, k2)
    
    kx, ky = unnormalize_kvector(graphene_lattice, [k1, k2, 0])
    t*sqrt(3+2*cos(sqrt(3)*kx*1.42)+4*cos(3/2*ky*1.42)*cos(sqrt(3)*kx/2*1.42))

end

function graphene_dos(t::T, mesh::R, histogram_width::S) where {T<:Number, R<:Number, S<:Number} 
    max_energy=3*abs(t)
    middle_index=round(Int, max_energy*histogram_width)+1
    num_indices=middle_index*2

    GrapheneDOS=zeros(num_indices)

    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]

    K=4*pi/(3*sqrt(3)*1.42);

    N=mesh
    for i in 1:N
        for j in 1:N
            kxnormal, kynormal=i/N, j/N
            kx, ky = unnormalize_kvector(graphene_lattice, [kxnormal, kynormal, 0])
            Ek=graphene_energy(t, kx, ky)


            GrapheneDOS[round(Int, histogram_width*Ek)+middle_index]=GrapheneDOS[round(Int, histogram_width*Ek)+middle_index]+(1/N)^2*histogram_width
            GrapheneDOS[-round(Int, histogram_width*Ek)+middle_index]=GrapheneDOS[-round(Int, histogram_width*Ek)+middle_index]+(1/N)^2*histogram_width
        end
    end
    return GrapheneDOS
end


function graphene_dos_quad(t::T, ϵ::R, δ::S; kwargs...) where {T<:Number, R<:Number, S<:Number}
    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]

    1/π*hcubature(vec->imag(-1/(ϵ-graphene_energy_normalizedk(t, graphene_lattice, vec[1], vec[2])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]
end


function check_graphene_dos_quad(t::T, δ::S, npoints::Int; kwargs...) where {T<:Number, S<:Number} 
    
    quaddos=[]
    for i in 1:npoints
        ω=i/npoints*abs(t)*3
        push!(quaddos, graphene_dos_quad(t, ω, δ; kwargs...))
    end
    
    return sum(quaddos*1/npoints*abs(t)*3)

end


"checks that the integrated value of the dos of graphene over all energies gives 2 orbitals per unit cell"
function check_graphene_dos(t::T, mesh::R, histogram_width::S) where {T<:Number, R<:Number, S<:Number} 

    graphene_dos_array = graphene_dos(t, mesh, histogram_width)

    return sum(graphene_dos_array*1/histogram_width)

end

function dirac_approximation_lower(k)
    return -6*k
end

function dirac_approximation_upper(k)
    return 6*k
end

function lower_band_integrand(k, theta, q, ω, delta)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
    
    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    return (2*mixedOverlap*(dirac_approximation_lower(k)-dirac_approximation_upper(kplusq))/((dirac_approximation_lower(k)-dirac_approximation_upper(kplusq))^2-(ω+1im*delta)^2))
end

function upper_band_integrand(k, theta, q, ω, delta)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    sameOverlap=1/2*(1+(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5 +delta/100000000))
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))

    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    a=2*sameOverlap*(dirac_approximation_upper(k)-dirac_approximation_upper(kplusq))/((dirac_approximation_upper(k)-dirac_approximation_upper(kplusq))^2-(ω+1im*delta)^2)
    b=2*mixedOverlap*(dirac_approximation_upper(k)-dirac_approximation_lower(kplusq))/((dirac_approximation_upper(k)-dirac_approximation_lower(kplusq))^2-(ω+1im*delta)^2)
    return (a+b)
end

function graphene_conductivity( μ, q, ω; kwargs... )
    delta=.01
    A= hcubature( x-> x[1]/(pi^2)*imag(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B=hcubature( x-> x[1]/(pi^2)*imag(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return -4im*ω/q^2*(B[1]+A[1])
end

function graphene_real_conductivity( μ, q, ω; kwargs... )
    delta=.01
    A= hcubature( x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B=hcubature( x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return 4*ω/q^2*(B[1]+A[1])
end

function graphene_epsilon( μ, q, ω; kwargs... )
    delta=.01
    A= hcubature( x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B=hcubature( x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return 1-90.5/q*(B[1]+A[1])
end

function find_graphene_plasmon(μ, q; nomegas=3, kwargs...)
    epsilon_array=Array{Float64, 1}(undef, nomegas)
    
    for i in 1:nomegas
        ω=i/nomegas*2μ
        epsilon_array[i]=(log∘abs)(graphene_epsilon( μ, q, ω; kwargs... ))
    end
    #return epsilon_array
    return argmin(epsilon_array)/nomegas*2μ
end



