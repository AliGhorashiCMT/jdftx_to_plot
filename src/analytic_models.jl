# We will include models for the dielectric function of graphene for future reference 

#We start first with the density of states of graphene


#=

The functions below are exact expressions for the dynamical polarization of graphene at charge neutrality 
and at finite doping. They are copied verbatim from: B Wunsch et al 2006 New J. Phys. 8 318


=#

function real_neutral(q::Real, w::Real)
    return -1*q^2/(4*(-w^2+36*q^2)^.5)
end

function imag_neutral(q::Real, w::Real)
    return -q^2/(4*(w^2-36*q^2)^.5)
end

function intraband_1a_real(q::Real, w::Real, μ::Real)
    return -2*μ/(36*pi)+1/4*(q^2)/sqrt(abs(36*(q^2)-w^2))
end

function gplus(x::Real)
    return x*sqrt(x^2-1)-log(x+sqrt(x^2-1))
end

function gminus(x::Real)
    if x>=0
        return x*sqrt(1-x^2)-atan(sqrt(1-x^2)/x)
    else
        return -pi+x*sqrt(1-x^2)+atan(sqrt(1-x^2)/(-x))
    end
end

function f(q::Real, w::Real)
    return 1/(4*pi)*q^2/sqrt(abs(w^2-(6*q)^2))
end

function intraband_1b_real(q::Real, w::Real, mu::Real)
    return -2*mu/(6^2*pi)+f(q, w)*(gplus((2*mu+w)/(6*q))-gplus((2*mu-w)/(6*q)))
end


function intraband_1b_imag(q::Real, w::Real, mu::Real)
    return f(q, w)*pi
end


function intraband_1a_imag(q::Real, w::Real, mu::Real)
    return f(q, w)*(gplus((2*mu-w)/(6*q))-gplus((2*mu+w)/(6*q)))
end


function intraband_2a_imag(q::Real, w::Real, mu::Real)
    return -f(q, w)*gplus((2*mu+w)/(6*q))
end


function intraband_2b_imag(q::Real, w::Real, mu::Real)
    return -f(q, w)*gminus((w-2*mu)/(6*q))
end

function intraband_2a_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)-f(q, w)*gminus((w-2*mu)/(6*q))
end

function intraband_2b_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)+f(q, w)*gplus((w+2*mu)/(6*q))
end

function intraband_3a_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)+f(q, w)*(gminus((2*mu+w)/(6*q))-gminus((w-2*mu)/(6*q)))
end

function intraband_3b_real(q::Real, w::Real, mu::Real)
    return -2*mu/(36*pi)+f(q, w)*(gplus((2*mu+w)/(6*q))-gplus((w-2*mu)/(6*q)))
end

function intraband_real_total(q::Real, w::Real, mu::Real)
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

function intraband_imag_total(q::Real, w::Real, mu::Real)
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

function graphene_total_polarization(q::Real, w::Real, mu::Real)

    return intraband_real_total(q, w, mu)+ (w<6q ? real_neutral(q, w) : 0 )

end

function graphene_total_impolarization(q::Real, w::Real, mu::Real)

    return intraband_imag_total(q, w, mu) + (w>6q ? imag_neutral(q, w) : 0 )

end


function exact_graphene_epsilon(q::Real, w::Real, mu::Real)
    return 1-e²ϵ/2/q*graphene_total_polarization(q, w, mu) 
end

function exact_graphene_plasmon(q::Real, mu::Real)
    Epsilons=zeros(1000)
    for i in 1:1000
        ω = mu*i/1000*3
        Epsilons[i] = 1-e²ϵ/2/q*graphene_total_polarization(q, ω, mu) 
    end
    return argmin(log.(abs.(Epsilons)))*3/1000*mu
end

function exact_graphene_plasmonq(ω::Real, mu::Real)
    logEpsilons=zeros(100000)
    for i in 1:100000
        q = mu*i/100000*3/6
        #logEpsilons[i] = log(abs(1-e²ϵ/2/q*(graphene_total_polarization(q, ω, mu))))
        logEpsilons[i] = log(abs(1-e²ϵ/2/q*(graphene_total_polarization(q, ω, mu)+graphene_total_impolarization(q, ω, mu))))
    end
    return argmin(logEpsilons)*3/100000*mu/6
end

function graphene_plasmon_confinement(λ::Real, μ::Real)
    ω=1.24/λ
    lambdaair=λ*1e-6
    lambdap=2*pi/exact_graphene_plasmonq(ω, μ)*1e-10
    return lambdaair/lambdap
end

function exact_graphene_landau_damping(q::Real, w::Real, δ::Real, mu::Real)
    RePolω = graphene_total_polarization(q, w, mu) 
    RePolδω = graphene_total_polarization(q, w+δ, mu)  
    ImPol = graphene_total_impolarization(q, w, mu)
    return ImPol*δ/(RePolδω-RePolω)
end

function exact_graphene_landau_damping(q::Real, δ::Real, mu::Real)
    plasmon = exact_graphene_plasmon(q, mu)
    RePolω = graphene_total_polarization(q, plasmon, mu) 
    RePolδω = graphene_total_polarization(q, plasmon+δ, mu)  
    ImPol = graphene_total_impolarization(q, plasmon, mu)
    return ImPol*δ/(RePolδω-RePolω)
end


function graphene_energy(t::Real, kx::Real, ky::Real)
    t*sqrt(3+2*cos(sqrt(3)*kx*1.42)+4*cos(3/2*ky*1.42)*cos(sqrt(3)*kx/2*1.42))
end

function graphene_energy_normalizedk(t::Real, graphene_lattice::Array{<:Array{<:Real, 1}, 1}, k1::Real, k2::Real)
    
    kx, ky = unnormalize_kvector(graphene_lattice, [k1, k2, 0])
    t*sqrt(3+2*cos(sqrt(3)*kx*1.42)+4*cos(3/2*ky*1.42)*cos(sqrt(3)*kx/2*1.42))

end

function graphene_dos(t::Real, mesh::Real, histogram_width::Real) 
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


function graphene_dos_quad(t::Real, ϵ::Real, δ::Real; kwargs...) 
    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]

    1/π*hcubature(vec->imag(-1/(ϵ-graphene_energy_normalizedk(t, graphene_lattice, vec[1], vec[2])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]
end


function check_graphene_dos_quad(t::Real, δ::Real, npoints::Int; kwargs...) 
    
    quaddos=[]
    for i in 1:npoints
        ω=i/npoints*abs(t)*3
        push!(quaddos, graphene_dos_quad(t, ω, δ; kwargs...))
    end
    
    return sum(quaddos*1/npoints*abs(t)*3)

end


"checks that the integrated value of the dos of graphene over all energies gives 2 orbitals per unit cell"
function check_graphene_dos(t::Real, mesh::Real, histogram_width::Real) 

    graphene_dos_array = graphene_dos(t, mesh, histogram_width)

    return sum(graphene_dos_array*1/histogram_width)

end

function dirac_approximation_lower(k::Real)
    return -6*k
end

function dirac_approximation_upper(k::Real)
    return 6*k
end

function lower_band_integrand(k::Real, theta::Real, q::Real, ω::Real, delta::Real)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))
    
    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    return (2*mixedOverlap*(dirac_approximation_lower(k)-dirac_approximation_upper(kplusq))/((dirac_approximation_lower(k)-dirac_approximation_upper(kplusq))^2-(ω+1im*delta)^2))
end

function upper_band_integrand(k::Real, theta::Real, q::Real, ω::Real, delta::Real)
    #Note that mixedOverlap has been changed to defy divide by zero error 
    sameOverlap=1/2*(1+(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5 +delta/100000000))
    mixedOverlap=1/2*(1-(k+q*cos(theta))/((k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000))

    kplusq=(k^2+q^2+2*k*q*cos(theta))^.5
    a=2*sameOverlap*(dirac_approximation_upper(k)-dirac_approximation_upper(kplusq))/((dirac_approximation_upper(k)-dirac_approximation_upper(kplusq))^2-(ω+1im*delta)^2)
    b=2*mixedOverlap*(dirac_approximation_upper(k)-dirac_approximation_lower(kplusq))/((dirac_approximation_upper(k)-dirac_approximation_lower(kplusq))^2-(ω+1im*delta)^2)
    return (a+b)
end

function graphene_conductivity( μ::Real, q::Real, ω::Real; kwargs... )
    delta=.01
    A= hcubature( x-> x[1]/(pi^2)*imag(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B=hcubature( x-> x[1]/(pi^2)*imag(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return -4im*ω/q^2*(B[1]+A[1])
end

function graphene_real_conductivity( μ::Real, q::Real, ω::Real; kwargs... )
    delta=.01
    A= hcubature( x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B=hcubature( x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return 4*ω/q^2*(B[1]+A[1])
end

function graphene_epsilon( μ::Real, q::Real, ω::Real; kwargs... )
    delta=.01
    A= hcubature( x-> x[1]/(pi^2)*real(lower_band_integrand(x[1], x[2], q, ω , delta)), [0, 0], [2, 2π]; kwargs...)
    B=hcubature( x-> x[1]/(pi^2)*real(upper_band_integrand(x[1], x[2], q, ω, delta)), [0, 0], [μ/6, 2π]; kwargs...)
    return 1-90.5/q*(B[1]+A[1])
end

function find_graphene_plasmon(μ::Real, q::Real; nomegas::Int=3, kwargs...)
    epsilon_array=Array{Float64, 1}(undef, nomegas)
    
    for i in 1:nomegas
        ω=i/nomegas*2μ
        epsilon_array[i]=(log∘abs)(graphene_epsilon( μ, q, ω; kwargs... ))
    end
    #return epsilon_array
    return argmin(epsilon_array)/nomegas*2μ
end


function graphene_electron_self_energy(ϵ::Real, μ::Real)
    abs(ϵ-μ)>0.2 ? 0.0183*abs(ϵ-sign(ϵ-μ)*0.2) : 0
end

function graphene_electron_real_self_energy(ϵ::Real, μ::Real)

    pyintegrate.quad(x-> -graphene_electron_self_energy(x, μ)/π, -8.4, 8.4,  wvar=ϵ, weight="cauchy", limit=1000, epsrel=1e-10, epsabs=1e-10)[1]

end

function graphene_analytic_real_self_energy(ϵ::Real, μ::Real)
    G=0.0183;
    w0=0.2;
    W=7;
    w=ϵ;
    return G/pi*(w0*log(real(abs((w+.001im+w0)^2 /(((w+.001im)-μ)^2-w0^2) ))) - w*log(real(abs(W^2*(w+.001im-μ+w0)/(((w+.001im)+w0)^2*(w+.001im-μ-w0))  )))  );

end

