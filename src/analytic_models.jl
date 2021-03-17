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

function exact_graphene_plasmon(q::Real, mu::Real; num_evals::Int= 1000, max_multiple_of_mu::Int=3)
    Epsilons=zeros(num_evals)
    for i in 1:num_evals
        ω = mu*i/num_evals*max_multiple_of_mu
        Epsilons[i] = 1-e²ϵ/2/q*graphene_total_polarization(q, ω, mu) 
    end
    return argmin(log.(abs.(Epsilons)))*max_multiple_of_mu/num_evals*mu
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


function marinko_graphene_landau_damping(q::Real, μ::Real; mesh::Int= 100, histrogram_width::Int=100)
    Marinko_Plasmon_Element=4π/137*6.6*3*100
    loss = 0
    plasmon = exact_graphene_plasmon(q, mu)
    for i in 1:mesh
        k=i/mesh*μ/2
        for j in 1:mesh
            theta=j/mesh*2*π
            kx, ky=k*cos(theta), k*sin(theta)
    
            kplusq=sqrt((kx+q)^2+ky^2)
    
            Eupperk, Elowerkplusq = dirac_approximation_upper(k), dirac_approximation_lower(kplusq)
            Eupperkplusq = dirac_approximation_upper(kplusq)

            delta=1
            overlapUL = 1-(k+q*cos(theta))/(Complex(k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000)
            overlapUL = 1/2*overlapUL;

            overlapUU = 1-(k+q*cos(theta))/(Complex(k^2+q^2+2*k*q*cos(theta))^.5+ delta/100000000)
            overlapUU = 1/2*overlapUU;

            fupperk = heaviside(μ-Eupperk)
            flowerkplusq = heaviside(μ-Elowerkplusq)
            fupperkplusq = heaviside(μ-Eupperkplusq)

            DiffEnergiesUL = Eupperk-Elowerkplusq
            DiffEnergiesUU = Eupperk-Eupperkplusq

            if abs(DiffEnergiesUL-plasmon)*histogram_width<0.5 && DiffEnergiesUL>0
                loss = loss + k*(flowerkplusq)*(1-fupperk)*overlapUL*Marinko_Plasmon_Element/q*plasmon*1/π^2*histogram_width*(μ/mesh*0.5)*(2π/mesh)
            end

            if abs(DiffEnergiesUU-plasmon)*histogram_width<0.5 && DiffEnergiesUU>0
                loss = loss + k*(fupperkplusq)*(1-fupperk)*overlapUU*Marinko_Plasmon_Element/q*plasmon*1/π^2*histrogram_width*(μ/mesh*0.5)*(2π/mesh)
            end
        end
    end
    return loss*2π/ħ
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

function graphene_dos_monte_carlo(t::Real, mesh::Real, histogram_width::Real) 
    max_energy=3*abs(t)
    middle_index=round(Int, max_energy*histogram_width)+1
    num_indices=middle_index*2

    GrapheneDOS=zeros(num_indices)

    a=1.42*sqrt(3)
    graphene_lattice=[[a, 0, 0], [-a/2, a*sqrt(3)/2, 0], [0, 0, 10]]

    K=4*pi/(3*sqrt(3)*1.42);

    N=mesh

    rand_kxs, rand_kys = rand(mesh), rand(mesh)

    for rkx in rand_kxs
        for rky in rand_kys
            kxnormal, kynormal = rkx, rky
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

"""
The matrix elements used in this function are taken from:

Park, Cheol-Hwan, et al. "Velocity renormalization and carrier lifetime in graphene from the electron-phonon interaction." Physical review letters 99.8 (2007): 086804.

"""
function graphene_numerical_self_energy(μ::Real; mesh1::Int=100, mesh2::Int=100, histogram_width::Real=100, NQs::Int=50)
    
    g = .035*13.605662285137 # The energy provided in the paper is given in Rydberg
    phononEnergy = 0.2 
    SelfEnergyMat=zeros(NQs)
        
    for ks in 1:NQs
        k=(ks-NQs/2)/NQs*0.4
        E=k*6
        print(ks); flush(stdout)
        
        for i in 1:mesh1
            for j in 1:mesh2
                q, theta=i/mesh1*1, j/mesh2*(2*π) 
                qx, qy=q*cos(theta), q*sin(theta)
                kplusq=sqrt((k+qx)^2+(qy)^2)
            
                for band in [1, 2]
                    
                    if band==1
                        Energy = dirac_approximation_upper(kplusq)
                    elseif band==2
                        Energy = dirac_approximation_lower(kplusq)
                    end
                    Occupation=heaviside(μ-Energy)
                    DiffEnergies1=Energy+phononEnergy
                    DiffEnergies2=Energy-phononEnergy
    
                    if  abs(E-DiffEnergies1)*histogram_width<.5
                        SelfEnergyMat[ks]=SelfEnergyMat[ks]+(1-Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                    end
                    if abs(E-DiffEnergies2)*histogram_width<.5
                        SelfEnergyMat[ks]=SelfEnergyMat[ks]+(Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                    end
                end
            end
        end
    end
    
    return SelfEnergyMat
end


function graphene_monte_carlo_self_energy(μ::Real; mesh1::Int=100, mesh2::Int=100, histogram_width::Real=100, NQs::Int=50)
    
    g = .035*13.605662285137 # The energy provided in the paper is given in Rydberg
    phononEnergy = 0.2 
    SelfEnergyMat=zeros(NQs)
        
    for ks in 1:NQs
        k=(ks-NQs/2)/NQs*0.4
        E=k*6
        print(ks); flush(stdout)
        
        random_ks = rand(mesh1)
        random_thetas = rand(mesh2)

        for rks in random_ks
            for thetas in random_thetas
                q, theta=1*rks, 2*π*thetas
                qx, qy=q*cos(theta), q*sin(theta)
                kplusq=sqrt((k+qx)^2+(qy)^2)
            
                for band in [1, 2]
                    
                    if band==1
                        Energy = dirac_approximation_upper(kplusq)
                    elseif band==2
                        Energy = dirac_approximation_lower(kplusq)
                    end
                    Occupation=heaviside(μ-Energy)
                    DiffEnergies1=Energy+phononEnergy
                    DiffEnergies2=Energy-phononEnergy
    
                    if  abs(E-DiffEnergies1)*histogram_width<.5
                        SelfEnergyMat[ks]=SelfEnergyMat[ks]+(1-Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                    end
                    if abs(E-DiffEnergies2)*histogram_width<.5
                        SelfEnergyMat[ks]=SelfEnergyMat[ks]+(Occupation)*q*π*g^2*(2*π/mesh1)*(1/mesh2)*histogram_width
                    end
                end
            end
        end
    end
    
    return SelfEnergyMat
end

function graphene_second_order_losses(; mesh1::Int=10, mesh2::Int=200, μ::Real=0.64, nlambda::Int=50, histogram_width::Real=100)
    OverallFactor=2π/ħ*(8π/137)*ħ^5*c^5/1e12*1/5.24
    DiffEnergies=zeros(nlambda)
    q=(μ/6)*(2/10)
    
    for lambdas in 1:nlambda
        
        lambda=3+lambdas/nlambda*6
        omega =1.24/lambda
    
        println("Plasmon Freq:", omega)
    
        flush(stdout)
        for i in 1:mesh1
            k=i/mesh1*.5
            for j in 1:mesh2
                theta=j/mesh2*2π
                kx, ky=k*cos(theta), k*sin(theta)
                kplusq=sqrt((kx+q)^2+ky^2)
                
                Ei = dirac_approximation_upper(k)
                Em = dirac_approximation_upper(kplusq)
                
                overlap=1
                f1=heaviside(μ-Ei)
                f2=1-heaviside(μ-Em)
                
                for i in 1:mesh1
                    k2=i/mesh1*.5
                    for j in 1:mesh2
                        theta2=j/mesh2*2π 
                        kplusqplusq2 = sqrt((kx+q+k2*cos(theta2))^2+(ky+k2*sin(theta2))^2)
                        kplusq2 = sqrt((kx+k2*cos(theta2))^2+(ky+k2*sin(theta2))^2)
                        Em2 = dirac_approximation_upper(kplusq2)

                        Ef=6*kplusqplusq2
                        f3=1-heaviside(μ-Ef)
                        
                        DiffEnergies2=Ef-Ei+0.2
                        if abs(DiffEnergies2-omega)*histogram_width<0.5 && DiffEnergies2>0
                            #k, k2 factors for area, f1, f2, f3 factors for occupation, k2 and omega factors for matrix element, the rest for integral 
                            DiffEnergies[lambdas]=DiffEnergies[lambdas]+k*k2*f3*f2*f1*overlap*abs(k/(Em-Ei-omega)+kplusq2/(Em2-Ei-0.2+0.01im))^2*1/omega*q*histogram_width*(1/mesh1*.5)^2*(2π/mesh2)^2*0.226576*OverallFactor
                        end
                    end
                end
            end
        end
    end
    return DiffEnergies
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

alevitov= 134/sqrt(3);
Klevitov=4*pi/(3*sqrt(3)*alevitov);
heaviside(x)= x>0 ? 1 : 0

#=
The functions below are to reproduce results from the following paper: 

Intrinsically undamped plasmon modes in narrow electron bands
Cyprian Lewandowski, Leonid Levitov
Proceedings of the National Academy of Sciences Oct 2019, 116 (42) 20869-20874; DOI: 10.1073/pnas.1909069116


=#

function levitov_epsilon(qx, qy, ω; kwargs...)
    q=sqrt(qx^2+qy^2)
    1-e²ϵ*1000/(2*q)*hcubature( x->levitov_integrand(x[1], x[1], qx, qy, ω, .1)*heaviside(limit_up_levitov(x[1])-x[2])*heaviside(-limit_dn_levitov(x[1])+x[2]), [-Klevitov, -Klevitov], [Klevitov, Klevitov]; kwargs...)[1]
end

function limit_dn_levitov(x)
    A1=heaviside(x+Klevitov)*heaviside(-x-Klevitov/2)*(-sqrt(3)*x-4*pi./(3*alevitov));
    A2=heaviside(x+Klevitov/2)*heaviside(-x+Klevitov/2)*(-4*pi/(6*alevitov));
    A3=heaviside(x-Klevitov/2)*heaviside(Klevitov-x)*(sqrt(3)*x-4*pi/(3*alevitov));
    A=A1+A2+A3;
end

function limit_up_levitov(x)
    B1=heaviside(x+Klevitov)*heaviside(-x-Klevitov/2)*(sqrt(3)*x+4*pi/(3*alevitov));
    B2=heaviside(x+Klevitov/2)*heaviside(-x+Klevitov/2)*(4*pi./(6*alevitov));
    B3=heaviside(x-Klevitov/2)*heaviside(Klevitov-x)*(-sqrt(3)*x+4*pi/(3*alevitov));
    B=B1+B2+B3;
end

function levitov_integrand(kx, ky, qx, qy, w, delta)
    kplusqy=ky+qy;
    kplusqx=kx+qx;
    E1=exp(alevitov*ky*1im)+exp(-(alevitov*kx*sqrt(3)/2+alevitov*ky/2)*1im)+exp((alevitov*kx*sqrt(3)/2-alevitov/2*ky)*1im);
    E2=exp(alevitov*kplusqy*1im)+exp(-(alevitov*kplusqx*sqrt(3)/2+alevitov*kplusqy/2)*1im)+exp((alevitov*kplusqx*sqrt(3)/2-alevitov/2*kplusqy)*1im);
    mixedOverlap=(1-cos(angle(E1)-angle(E2)))/2;
    sameOverlap=(1+cos(angle(E1)-angle(E2)))/2;
    Up1=3.75/3*abs(E1);
    Up2=3.75/3*abs(E2);
    a=2*sameOverlap*(Up1-Up2)/((Up1-Up2)^2-(w+1im*delta)^2);
    b=2*mixedOverlap*(Up1+Up2)/((Up1+Up2)^2-(w+1im*delta)^2);
    fullbands1= 1/(12.12*pi^2)*(a+b)*heaviside(1.81-Up1);
    fullbands2 = 2/(12.12*pi^2)*mixedOverlap*(-Up1-Up2)/((Up1+Up2)^2-(w+1im*delta)^2);
    fullbands=fullbands1+fullbands2;
    return fullbands
end

